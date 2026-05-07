#!/usr/bin/env Rscript

#This script export "unassigned" cells from Seurat object (SeuratV5)
#into MTX + genes.txt+ cells,txt to be loaded in Python and annotated with cell typist 


suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(future)
})

start_time <- Sys.time()

# Define the arguments for bash 
option_list <- list(
  make_option(c("--seurat_rds"), type = "character", help = "Path to Seurat .rds", metavar = "FILE"),
  make_option(c("--out_dir"),    type = "character", default = "unassigned_only_final",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("--assay"),      type = "character", default = "RNA",
              help = "Assay to use if present (fallback to DefaultAssay if not) [default: %default]"),
  make_option(c("--layer_prefix"), type = "character", default = "data",
              help = "Layer family to join (e.g. 'data' joins data.1,data.2,...) [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer", default = 4,
              help = "Parallel workers (for future plan) [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate the input (rds.seurat)
stopifnot(!is.null(opt$seurat_rds))
stopifnot(file.exists(opt$seurat_rds))

#Set seed and workers
set.seed(opt$seed) 
plan(multicore, workers = opt$workers)

# Output dir
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("# Seurat rds: ", opt$seurat_rds, "\n")
cat("# Out dir:   ", normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("# Assay pref:", opt$assay, "\n")
cat("# Layer pref:", opt$layer_prefix, "\n")
cat("# Seed:      ", opt$seed, "\n")
cat("# Workers:   ", opt$workers, "\n\n")

#  Load object 
cat("# Loading Seurat object...\n")
obj <- readRDS(opt$seurat_rds)
cat("# Cells in Seurat: ", ncol(obj), "\n")

# 1:Identify unassigned cells (strictly by cell_type == NA)

meta_cols <- colnames(obj@meta.data)

if (!("cell_type" %in% meta_cols)) {
  stop("Required column 'cell_type' not found in meta.data")
}

un_cells <- colnames(obj)[is.na(obj$cell_type)]

cat("# Unassigned cells:", length(un_cells), "\n")
if (length(un_cells) == 0) {
  stop("No unassigned cells found. Nothing to export.")
}


#Only use the unassigned cells 

obj_un <- subset(obj, cells = un_cells)

#  2 : Chosse the assay to export , in this case are layers  
assay_use <- if (opt$assay %in% Assays(obj_un)) opt$assay else DefaultAssay(obj_un)
DefaultAssay(obj_un) <- assay_use

cat("# Using assay: ", assay_use, "\n")
cat("# Layers BEFORE join:\n")
print(Layers(obj_un[[assay_use]]))

# 3: Join layers (e.g., data.* -> data)  -> Because of the integration we want to merge all the layers
cat("# Joining layers: ", opt$layer_prefix, ".* -> ", opt$layer_prefix, "\n", sep = "")
obj_un[[assay_use]] <- JoinLayers(obj_un[[assay_use]], layers = opt$layer_prefix)

cat("# Layers AFTER join:\n")
print(Layers(obj_un[[assay_use]]))

# 4: Extract the merged layer matrix and export to MTX+ gene/cell.txt
cat("# Extracting layer: ", opt$layer_prefix, "\n", sep = "")
mat <- LayerData(obj_un[[assay_use]], layer = opt$layer_prefix)

cat("# Matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "), "\n", sep = "")

mtx_path   <- file.path(opt$out_dir, "unassigned_logNorm.mtx")
genes_path <- file.path(opt$out_dir, "genes.txt")
cells_path <- file.path(opt$out_dir, "cells.txt")

#Save the matrixes 
Matrix::writeMM(mat, mtx_path)
writeLines(rownames(mat), genes_path)
writeLines(colnames(mat), cells_path)

cat("\n# Wrote:\n")
cat(" - ", normalizePath(mtx_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(" - ", normalizePath(genes_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(" - ", normalizePath(cells_path, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("# Files exist? ",
    file.exists(mtx_path), file.exists(genes_path), file.exists(cells_path), "\n")

end_time <- Sys.time()
cat("\n# DONE. Time taken: ", end_time - start_time, "\n")
