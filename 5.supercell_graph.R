#!/usr/bin/env Rscript

# Super cell adjustments
# Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(optparse)
})

# Arguments
option_list <- list(
  make_option(c("-s", "--seurat_rds"),
              type    = "character",
              help    = "Path to annotated Seurat object (.rds)"),
  make_option(c("-o", "--out_dir"),
              type    = "character",
              default = "supercell_results",
              help    = "Output directory [default: %default]"),
  make_option(c("-g", "--gammas"),
              type    = "character",
              default = "10,20,30,50",
              help    = "Graining levels comma-separated [default: %default]"),
  make_option(c("--cell_type_col"),
              type    = "character",
              default = "cell_type",
              help    = "Column name for cell type annotation [default: %default]"),
  make_option(c("--phenotype_col"),
              type    = "character",
              default = "is_AD",
              help    = "Column name for phenotype [default: %default]"),
  make_option(c("--reduction"),
              type    = "character",
              default = "harmony",
              help    = "Reduction to use [default: %default]"),
  make_option(c("--n_pc"),
              type    = "integer",
              default = 30,
              help    = "Number of PCs [default: %default]"),
  make_option(c("--k_knn"),
              type    = "integer",
              default = 5,
              help    = "Number of nearest neighbors [default: %default]"),
  make_option(c("--do_approx"),
              action  = "store_true",
              default = FALSE,
              help    = "Use approximate kNN [default: %default]"),
  make_option(c("--approx_n"),
              type    = "integer",
              default = 20000,
              help    = "Cells for approximate kNN [default: %default]"),
  make_option(c("--seed"),
              type    = "integer",
              default = 42,
              help    = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validation
stopifnot(file.exists(opt$seurat_rds))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("CONSTRUCCIÓN DE METACÉLULAS\n")
cat("Objeto:        ", opt$seurat_rds, "\n")
cat("Salida:        ", opt$out_dir, "\n")
cat("Reducción:     ", opt$reduction, "\n")
cat("n_pc:          ", opt$n_pc, "\n")
cat("k_knn:         ", opt$k_knn, "\n")
cat("cell_type_col: ", opt$cell_type_col, "\n")
cat("phenotype_col: ", opt$phenotype_col, "\n")
cat("Gammas:        ", opt$gammas, "\n")
cat("Inicio:        ", format(Sys.time()), "\n\n")

# Load Seurat
cat("[1/4] Cargando Seurat...\n")
obj <- readRDS(opt$seurat_rds)
cat(sprintf("      %s células | %s genes\n",
            format(ncol(obj), big.mark = ","),
            format(nrow(obj), big.mark = ",")))

# Filter NAs
cat("[2/4] Filtrando células sin anotación...\n")
cells_keep <- !is.na(obj@meta.data[[opt$cell_type_col]]) &
  !is.na(obj@meta.data[[opt$phenotype_col]]) &
  obj@meta.data[[opt$cell_type_col]] != "" &
  obj@meta.data[[opt$phenotype_col]] != ""

obj_clean <- obj[, cells_keep]
cat(sprintf("      Células limpias: %s\n", format(sum(cells_keep), big.mark = ",")))

# Extract embedding — UNA sola vez, fuera del loop
cat("[3/4] Extrayendo embedding...\n")
embedding <- Embeddings(obj_clean, reduction = opt$reduction)[, 1:opt$n_pc]
cat(sprintf("      %s células × %s dimensiones\n",
            nrow(embedding), ncol(embedding)))

# Parse gammas
gammas <- as.numeric(strsplit(opt$gammas, ",")[[1]])
cat(sprintf("\n[4/4] Corriendo SCimplify_from_embedding para %s gammas...\n",
            length(gammas)))

# Loop por gamma — SCimplify corre una vez por cada gamma
for (gamma in gammas) {
  
  cat(sprintf("\n  [gamma=%s] Construyendo metacélulas...\n", gamma))
  
  SC <- SCimplify_from_embedding(
    X                    = embedding,
    cell.annotation      = obj_clean@meta.data[[opt$cell_type_col]],
    cell.split.condition = obj_clean@meta.data[[opt$phenotype_col]],
    gamma                = gamma,
    k.knn                = opt$k_knn,
    n.pc                 = opt$n_pc,
    do.approx            = opt$do_approx,
    approx.N             = opt$approx_n,
    seed                 = opt$seed
  )
  
  cat(sprintf("  [gamma=%s] Metacélulas construidas: %s\n",
              gamma, format(SC$N.SC, big.mark = ",")))
  cat(sprintf("  [gamma=%s] Células asignadas:       %s\n",
              gamma, length(SC$membership)))
  
  # Un .rds por gamma
  out_rds <- file.path(opt$out_dir, sprintf("SC_gamma%s.rds", gamma))
  saveRDS(SC, out_rds)
  cat(sprintf("  [gamma=%s] Guardado: %s\n", gamma, out_rds))
}

cat(sprintf("\n✓ Fin: %s\n", format(Sys.time())))

