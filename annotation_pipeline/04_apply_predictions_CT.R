#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(future)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(RColorBrewer)
})

start_time <- Sys.time()

#arguments
option_list <- list(
  make_option(c("--seurat_rds"), type = "character", help = "Input Seurat .rds", metavar = "FILE"),
  make_option(c("--pred_csv"),   type = "character", help = "CellTypist output CSV (unassigned_celltypist_predictions.csv)", metavar = "FILE"),
  make_option(c("--out_dir"),    type = "character", default = "annotation_with_celltypist",
              help = "Output folder for annotated RDS + PNGs [default: %default]", metavar = "DIR"),
  make_option(c("--out_rds"),    type = "character", default = "merged_harmony_integrated_annotated_plus_celltypist.rds",
              help = "Output RDS filename (saved inside out_dir) [default: %default]", metavar = "FILE"),
  make_option(c("--min_prob"),   type = "double", default = NA,
              help = "Optional: if provided, only fill cell_type when celltypist_max_prob >= min_prob", metavar = "FLOAT"),
  make_option(c("-w", "--workers"), type = "integer", default = 4,
              help = "Parallel workers (for future plan) [default: %default]"),
  make_option(c("--seed"),       type = "integer", default = 42,
              help = "Random seed (mostly for reproducibility of plotting) [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
#Validate the existance of the rds and csv
stopifnot(!is.null(opt$seurat_rds), file.exists(opt$seurat_rds))
stopifnot(!is.null(opt$pred_csv),   file.exists(opt$pred_csv))

#workers and seed
plan(multicore, workers = opt$workers)
set.seed(opt$seed)

#create the carpet , folder out 
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("# Seurat rds: ", opt$seurat_rds, "\n")
cat("# Pred CSV:   ", opt$pred_csv, "\n")
cat("# Out dir:    ", normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("# Out rds:    ", opt$out_rds, "\n")
cat("# min_prob:   ", ifelse(is.na(opt$min_prob), "NA (no extra filter)", opt$min_prob), "\n")
cat("# Workers:   ", opt$workers, "\n\n")
cat("# Seed:       ", opt$seed, "\n\n")

# Inputs 
cat("# Loading Seurat object...\n")
obj <- readRDS(opt$seurat_rds)
cat("# Cells in Seurat: ", ncol(obj), "\n")

#Load the csv with the predictions of cel typist 
cat("# Loading CellTypist predictions...\n")
#pred: the data frame that came from Cell Typist 
pred <- readr::read_csv(opt$pred_csv, show_col_types = FALSE) 

# Ensure types matches  , normalize the input 
# Ensure key columns are stored as character rather than factor.
pred <- pred %>% #apply in "pred"
  mutate(
    cell = as.character(cell),
    celltypist_label = as.character(celltypist_label),
    celltype_simple = as.character(celltype_simple)
  )

cat("# Pred rows: ", nrow(pred), "\n")
cat("# Unique cells in pred: ", length(unique(pred$cell)), "\n")

# We need to match the cells 
# The merge , here is the importante of the match beetwen cells in the original Seurat object and the predictions
# colnames(obj) : Seurat , pred$cell : CellTypist
common_cells <- intersect(colnames(obj), pred$cell)
cat("# Cells that match Seurat: ", length(common_cells), "\n") #Cells that match 
cat("# Cells in pred not in Seurat: ", sum(!pred$cell %in% colnames(obj)), "\n") #Cells that are in prediction.csv but not in seurat , ideal =0
cat("# Cells in Seurat not in pred: ", sum(!colnames(obj) %in% pred$cell), "\n") # Cells only in Seurat -> assigned the first time 

if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat object and prediction CSV.")
}

# Create this new colums in Seurat in meta.data  #is the same as ->  obj@meta.data$nombre_columna
#Create one row per cell
# NA :text
obj$celltypist_label <- NA_character_  
obj$celltypist_simple <- NA_character_
obj$celltypist_max_prob <- NA_real_
obj$celltypist_margin <- NA_real_
obj$assigned_by_celltypist <- FALSE

# map -> Como un vector nombrado , tipo diccionario 
# Diferents label are associated to a barcode 

lab_map <- setNames(pred$celltypist_label, pred$cell)  #specific type
sim_map <- setNames(pred$celltype_simple, pred$cell) #symple type 
prob_map <- setNames(pred$max_prob, pred$cell) # prob
mar_map  <- setNames(pred$margin_top1_top2, pred$cell) #margen

# coomon cells : is a vector of barcodes 
#lab_map -> give the name of each label 
obj$celltypist_label[common_cells]  <- lab_map[common_cells]
obj$celltypist_simple[common_cells] <- sim_map[common_cells]
obj$celltypist_max_prob[common_cells] <- prob_map[common_cells]
obj$celltypist_margin[common_cells]   <- mar_map[common_cells]

# assigned_by_celltypist: label exists and is not "Unassigned"
#TRUE -> Cell typist could assigned 
#False -> Cell typist could not assigned 
obj$assigned_by_celltypist[common_cells] <- !is.na(obj$celltypist_simple[common_cells]) &
  obj$celltypist_simple[common_cells] != "Unassigned"

# Optional additional confidence filter before filling cell_type
#You have to put the min conf of probability , if not is only select the cell wich dont have "NA"
if (!is.na(opt$min_prob)) {
  obj$assigned_by_celltypist[common_cells] <- obj$assigned_by_celltypist[common_cells] &
    !is.na(obj$celltypist_max_prob[common_cells]) &
    obj$celltypist_max_prob[common_cells] >= opt$min_prob
}

cat("# assigned_by_celltypist table:\n")
print(table(obj$assigned_by_celltypist))

#  Fill NA cell_type only 
if (!("cell_type" %in% colnames(obj@meta.data))) {
  stop("Required column 'cell_type' not found in Seurat meta.data")
}
#What cell did not have a previos annotation ? select the cell NA in the original anotation
#True: cell_type was NA / False : cell_type with previos anno
#na: cells NA in the original anno
na_cells <- colnames(obj)[is.na(obj$cell_type)]
cat("# Original NA cells: ", length(na_cells), "\n")

#Intersect common cells who was na_cells and CT annotate
fill_cells <- intersect(na_cells, common_cells)
#Make the vector TRUE/FALSE , just keep true (not previos assigned a cell type)
fill_cells <- fill_cells[obj$assigned_by_celltypist[fill_cells]]  # only those confidently assigned

#fill the cell_type just in that cells 
obj$cell_type[fill_cells] <- obj$celltypist_simple[fill_cells]
#for plots 
obj$cell_type_plot_final <- ifelse(is.na(obj$cell_type), "Unassigned", obj$cell_type)

cat("# Filled: ", length(fill_cells), "\n") #How much could be filled with CT?
cat("# Remaining NA: ", sum(is.na(obj$cell_type)), "\n") # Remain NA?

# Track annotation source, where this cell was annoateted?
obj$annotation_source <- "original"
obj$annotation_source[fill_cells] <- "celltypist"
obj$annotation_source[colnames(obj)[is.na(obj$cell_type)]] <- "unassigned"

cat("# annotation_source table:\n")
print(table(obj$annotation_source))

#  Plots :) 

# Note: ggrepel overlap control is global option (DimPlot doesn't accept max.overlaps arg directly)
options(ggrepel.max.overlaps = 500)

# helper to save plots
save_png <- function(p, filename, width = 2400, height = 2000, res = 300) {
  png(file.path(opt$out_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# 1: Full UMAP: cell_type_plot_final 
p1 <- DimPlot(obj, group.by = "cell_type_plot_final", label = TRUE) + NoLegend()
save_png(p1, "umap_cell_type_plot_final_basic.png")

# 2: Full UMAP
p2 <- DimPlot(obj,
              group.by = "cell_type_plot_final",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.25,
              label.size = 3.5) + NoLegend()
save_png(p2, "umap_cell_type_plot_final_repel_pt025_label35.png")

p3 <- DimPlot(obj,
              group.by = "cell_type_plot_final",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.2)
save_png(p3, "umap_cell_type_plot_final_repel_pt02.png")

# 3: annotation_source
p4 <- DimPlot(obj, group.by = "annotation_source")
save_png(p4, "umap_annotation_source.png")

# 4:FeaturePlot: max_prob
p5 <- FeaturePlot(obj, features = "celltypist_max_prob")
save_png(p5, "featureplot_celltypist_max_prob.png")

# UMAP just with 6 major cell types 
types_keep <- c(
  "Excitatory Neurons",
  "Inhibitory Neurons",
  "Astrocyte",
  "OPCs",
  "Oligodendrocytes",
  "Microglia"
)

obj_sub <- subset(obj, subset = cell_type_plot_final %in% types_keep)

# enforce factor order (important for consistent coloring)
obj_sub$cell_type_plot_final <- factor(obj_sub$cell_type_plot_final, levels = types_keep)

p6 <- DimPlot(obj_sub,
              group.by = "cell_type_plot_final",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.25,
              label.size = 3.5) + NoLegend()
save_png(p6, "umap_types_keep_only_basic.png")

# 
p7 <- DimPlot(obj_sub,
              group.by = "cell_type_plot_final",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.2,
              label.size = 3) +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  ggtitle("Major Brain Cell Types") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  NoLegend()
save_png(p7, "umap_types_keep_only_brewer_Set2.png")

#  custom colors attempt
custom_colors <- c(
  "Excitatory Neurons" = "#6495ED",
  "Inhibitory Neurons" = "#FF6347",
  "Astrocyte"          = "#9ACD32",
  "OPCs"               = "#9370DB",
  "Oligodendrocytes"   = "#FFA500",
  "Microglia"          = "#FFD700"
)

p8 <- DimPlot(obj_sub,
              group.by = "cell_type_plot_final",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.25,
              label.size = 3.5) +
  scale_color_manual(values = custom_colors) +
  ggtitle("Major Brain Cell Types") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  NoLegend()
save_png(p8, "umap_types_keep_only_custom_colors.png")

#  Keep-only as NA for others (white background) 
obj$cell_type_keep_only <- ifelse(obj$cell_type_plot_final %in% types_keep,
                                  obj$cell_type_plot_final,
                                  NA)

p9 <- DimPlot(obj,
              group.by = "cell_type_keep_only",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.25,
              label.size = 3.5,
              na.value = "white") + NoLegend()
save_png(p9, "umap_keep_only_na_white.png")

# Save final Seurat object
out_rds_path <- file.path(opt$out_dir, opt$out_rds)
saveRDS(obj, out_rds_path)
cat("\n# Saved annotated Seurat RDS: ", normalizePath(out_rds_path, winslash = "/", mustWork = FALSE), "\n")

end_time <- Sys.time()
cat("# DONE. Time taken: ", end_time - start_time, "\n")
