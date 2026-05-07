start_time <- Sys.time()

#Define option list for inputs 
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(future)
})

#Define option list for inputs 
option_list <- list(
  make_option(c("--seurat_rds"), type = "character", help = "Ruta al Seurat integrado (.rds)", metavar = "FILE"),
  make_option(c("--anno_csv"),   type = "character", help = "Ruta al CSV de anotación celular", metavar = "FILE"),
  make_option(c("--out_dir"),    type = "character", default = "out_annotated_final", help = "Directorio de salida [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer",default = 4 ,help = "Parallel workers [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,help = "Random seed [default: %default]")
  )

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

stopifnot(!is.null(opt$seurat_rds), !is.null(opt$anno_csv))
stopifnot(file.exists(opt$seurat_rds), file.exists(opt$anno_csv))

# Parallelization

plan(multicore, workers = opt$workers)
options(future.globals.maxSize = 200 * 1024^3)


# Output dir

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(opt$out_dir)

cat("# Seurat rds: ", opt$seurat_rds, "\n")
cat("# Annotation csv: ", opt$anno_csv, "\n")
cat("# Out dir: ", normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE), "\n")

# Cargar inputs

cat("# Loading annotation CSV...\n")
celular_annotation <- readr::read_csv(opt$anno_csv, show_col_types = FALSE)

cat("# Loading Seurat object...\n")
seurat_list_integrated <- readRDS(opt$seurat_rds)

cat("# Cells in Seurat: ", ncol(seurat_list_integrated), "\n")


# Limpiar annotation (solo columnas clave)

cat("# Cleaning annotation...\n")
celular_annotation_clean <- celular_annotation %>%
  dplyr::select(cell, cell.type, state) %>%
  dplyr::distinct()

annotation_map <- celular_annotation_clean %>%
  tibble::column_to_rownames("cell")


# Crear llave cell_key para empatar -> En Seurat: ACTACGA...-1_1 y En CSV:   190403-B4-A_ACTACGA...-1
# Quitamos _1 y agregamos batch al inicio
cat("# Building cell_key...\n")

barcode_clean <- sub("_[0-9]+$", "", colnames(seurat_list_integrated))  # quita _1, _2...
batch <- seurat_list_integrated$libraryBatch                            # batch por célula
cell_key <- paste0(batch, "_", barcode_clean)


# Pegar anotación

cat("# Annotating Seurat metadata...\n")
seurat_list_integrated$cell_type  <- annotation_map[cell_key, "cell.type", drop = TRUE]
seurat_list_integrated$cell_state <- annotation_map[cell_key, "state", drop = TRUE]

# Para graficar NA como "Unassigned"
seurat_list_integrated$cell_type_plot <- ifelse(
  is.na(seurat_list_integrated$cell_type),
  "Unassigned",
  seurat_list_integrated$cell_type
)


# Resumen de anotación

total_cells <- ncol(seurat_list_integrated)
annotated_cells <- sum(!is.na(seurat_list_integrated$cell_type))
unannotated_cells <- sum(is.na(seurat_list_integrated$cell_type))

cat("\n# ===== Annotation summary =====\n")
cat("Total cells: ", total_cells, "\n")
cat("Annotated cells: ", annotated_cells, "\n")
cat("Unassigned cells (NA): ", unannotated_cells, "\n")
cat("Percent annotated: ", round(annotated_cells / total_cells * 100, 2), "%\n")
cat("Percent unassigned: ", round(unannotated_cells / total_cells * 100, 2), "%\n")

# Guardar un txt de resumen
writeLines(
  c(
    paste0("Total cells: ", total_cells),
    paste0("Annotated cells: ", annotated_cells),
    paste0("Unassigned cells (NA): ", unannotated_cells),
    paste0("Percent annotated: ", round(annotated_cells / total_cells * 100, 2), "%"),
    paste0("Percent unassigned: ", round(unannotated_cells / total_cells * 100, 2), "%")
  ),
  con = "annotation_summary.txt"
)

# Diagnóstico: ejemplos de NA y batches
cat("\n# Example NA keys:\n")
print(head(cell_key[is.na(seurat_list_integrated$cell_type)]))

cat("\n# NA per libraryBatch:\n")
print(table(seurat_list_integrated$libraryBatch[is.na(seurat_list_integrated$cell_type)]))


# Exportar outputs

cat("\n# Saving outputs...\n")

# Metadata completa anotada
write.csv(
  seurat_list_integrated@meta.data,
  file = "merged_harmony_integrated_cell_metadata_annotated.csv",
  row.names = TRUE
)

# Seurat anotado
saveRDS(
  seurat_list_integrated,
  file = "merged_harmony_integrated_annotated.rds"
)

# UMAPs anotados

cat("# Writing UMAP PDFs...\n")

pdf("umap_cell_type.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type", label = TRUE) + NoLegend())
dev.off()

pdf("umap_cell_type_with_unassigned.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type_plot", label = TRUE) + NoLegend())
dev.off()

pdf("umap_cell_state.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_state", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

end_time <- Sys.time()
cat("\n# DONE. Time taken: ", end_time - start_time, "\n")
