#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(optparse)
  library(future)
})

# Options
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
              default = "10,20,50",
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
              help    = "Reduction to use for embedding [default: %default]"),
  make_option(c("--n_pc"),
              type    = "integer",
              default = 30,
              help    = "Number of PCs to use [default: %default]"),
  make_option(c("--k_knn"),
              type    = "integer",
              default = 5,
              help    = "Number of nearest neighbors [default: %default]"),
  make_option(c("--do_approx"),
              action  = "store_true",
              default = FALSE,
              help    = "Use approximate kNN for large datasets [default: %default]"),
  make_option(c("--approx_n"),
              type    = "integer",
              default = 100000,
              help    = "Number of cells for approximate kNN [default: %default]"),
  make_option(c("-w", "--workers"),
              type    = "integer",
              default = 4,
              help    = "Parallel workers [default: %default]"),
  make_option(c("--seed"),
              type    = "integer",
              default = 42,
              help    = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# validacion
stopifnot(file.exists(opt$seurat_rds))

# Parsear gammas
gammas <- as.numeric(strsplit(opt$gammas, ",")[[1]])
cat("Gammas a correr:", paste(gammas, collapse=", "), "\n")

# Crear directorio de salida
fecha   <- format(Sys.Date(), "%Y-%m-%d")
out_dir <- opt$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Paralelización
plan(multicore, workers = opt$workers)
options(future.globals.maxSize = 200 * 1024^3)

cat("CONSTRUCCIÓN DE METACÉLULAS CON SUPERCELL\n")
cat("Objeto:       ", opt$seurat_rds, "\n")
cat("Salida:       ", out_dir, "\n")
cat("Gammas:       ", paste(gammas, collapse=", "), "\n")
cat("Reducción:    ", opt$reduction, "\n")
cat("Cell type col:", opt$cell_type_col, "\n")
cat("Phenotype col:", opt$phenotype_col, "\n")
cat("Workers:      ", opt$workers, "\n")
cat("Inicio:       ", format(Sys.time()), "\n")


#  1. CARGAR OBJ
cat("[1/4] Cargando objeto Seurat...\n")
obj <- readRDS(opt$seurat_rds)
cat(sprintf("      OK %s células | %s genes\n",
            format(ncol(obj), big.mark=","),
            format(nrow(obj), big.mark=",")))

# 2. FILTRAR NAs 
cat("[2/4] Filtrando células sin anotación...\n")

cells_keep <- !is.na(obj@meta.data[[opt$cell_type_col]]) &
  !is.na(obj@meta.data[[opt$phenotype_col]]) &
  obj@meta.data[[opt$cell_type_col]] != "" &
  obj@meta.data[[opt$phenotype_col]] != ""

cat(sprintf("      Células originales:       %s\n", format(ncol(obj), big.mark=",")))
cat(sprintf("      Células filtradas (NA):   %s\n", format(sum(!cells_keep), big.mark=",")))
cat(sprintf("      Células para SuperCell:   %s\n", format(sum(cells_keep), big.mark=",")))

obj_clean <- obj[, cells_keep]

# Extraer embedding y anotaciones
cat(sprintf("[3/4] Extrayendo embedding '%s'...\n", opt$reduction))
harmony_embedding <- Embeddings(obj_clean, reduction = opt$reduction)
cell_type <- droplevels(as.factor(obj_clean@meta.data[[opt$cell_type_col]]))
fenotipo  <- droplevels(as.factor(obj_clean@meta.data[[opt$phenotype_col]]))

cat(sprintf("      Tipos celulares: %s\n", paste(levels(cell_type), collapse=", ")))
cat(sprintf("      Fenotipos:       %s\n", paste(levels(fenotipo),  collapse=", ")))

#  3. LOOP DE GAMMAS - Granulacion  
cat("\n[4/4] Construyendo metacélulas...\n")

resumen <- data.frame()

for(gamma in gammas){
  
  cat(sprintf("\n  ── Gamma = %d \n", gamma))
  
  SC <- SCimplify_from_embedding(
    X                    = harmony_embedding,
    gamma                = gamma,
    k.knn                = opt$k_knn,
    n.pc                 = opt$n_pc,
    do.approx            = opt$do_approx,
    approx.N             = opt$approx_n,
    seed                 = opt$seed,
    cell.annotation      = cell_type,
    cell.split.condition = fenotipo
  )
  
  cat(sprintf("  Metacélulas construidas: %d\n", SC$N.SC))
  cat(sprintf("  Tamaño promedio:         %.1f células/metacélula\n",
              mean(SC$supercell_size)))
  
  # Guardar objeto
  nombre_rds <- file.path(out_dir,
                          paste0("SC_gamma", gamma, "_", fecha, ".rds"))
  saveRDS(SC, file = nombre_rds)
  cat(sprintf("  Guardado: %s\n", nombre_rds))
  
  # Plot por tipo celular
  png(file.path(out_dir, paste0("plot_celltype_gamma", gamma, "_", fecha, ".png")),
      width=1400, height=1000, res=150)
  supercell_plot(SC$graph.supercells,
                 group = SC$SC.cell.annotation.,
                 main  = paste("Cell type | gamma =", gamma),
                 seed  = opt$seed)
  dev.off()
  
  # Plot por fenotipo
  png(file.path(out_dir, paste0("plot_phenotype_gamma", gamma, "_", fecha, ".png")),
      width=1400, height=1000, res=150)
  supercell_plot(SC$graph.supercells,
                 group = SC$SC.cell.split.condition.,
                 main  = paste("Phenotype | gamma =", gamma),
                 seed  = opt$seed)
  dev.off()
  
  # Tabla resumen de metacélulas por tipo celular y fenotipo
  tab <- as.data.frame(table(
    cell_type = SC$SC.cell.annotation.,
    phenotype = SC$SC.cell.split.condition.
  ))
  write.csv(tab,
            file.path(out_dir, paste0("summary_gamma", gamma, "_", fecha, ".csv")),
            row.names = FALSE)
  
  # Acumular resumen general
  resumen <- rbind(resumen, data.frame(
    gamma            = gamma,
    n_metacells      = SC$N.SC,
    mean_size        = round(mean(SC$supercell_size), 1),
    min_size         = min(SC$supercell_size),
    max_size         = max(SC$supercell_size)
  ))
}

# Guardar resumen general
write.csv(resumen,
          file.path(out_dir, paste0("resumen_gammas_", fecha, ".csv")),
          row.names = FALSE)


cat("RESUMEN FINAL\n")
print(resumen)
cat(sprintf("\nArchivos guardados en: %s\n", out_dir))
cat(sprintf("Fin: %s\n", format(Sys.time())))