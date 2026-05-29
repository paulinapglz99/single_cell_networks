#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Matrix)
  library(optparse)
})

# ── OPCIONES ─────────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-s", "--seurat_rds"),
              type    = "character",
              help    = "Path to annotated Seurat object (.rds)"),
  make_option(c("-c", "--sc_rds"),
              type    = "character",
              help    = "Path to SuperCell object (.rds)"),
  make_option(c("-o", "--out_dir"),
              type    = "character",
              default = "matrices_aracne",
              help    = "Output directory [default: %default]"),
  make_option(c("--cell_type_col"),
              type    = "character",
              default = "cell_type",
              help    = "Column name for cell type annotation [default: %default]"),
  make_option(c("--phenotype_col"),
              type    = "character",
              default = "is_AD",
              help    = "Column name for phenotype [default: %default]"),
  make_option(c("--cell_types"),
              type    = "character",
              default = "Excitatory Neurons,Inhibitory Neurons,Astrocyte,Oligodendrocytes,Endothelial,OPCs,Microglia",
              help    = "Comma-separated cell types to keep [default: %default]"),
  make_option(c("--slot"),
              type    = "character",
              default = "data",
              help    = "Seurat slot to use: counts or data [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── VALIDACION ────────────────────────────────────────────────────────────────
stopifnot(file.exists(opt$seurat_rds))
stopifnot(file.exists(opt$sc_rds))

# ── SETUP ─────────────────────────────────────────────────────────────────────
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cell_types_keep <- trimws(strsplit(opt$cell_types, ",")[[1]])

cat("BUILD MATRICES PARA ARACNE\n")
cat("Seurat:       ", opt$seurat_rds, "\n")
cat("SC object:    ", opt$sc_rds, "\n")
cat("Salida:       ", opt$out_dir, "\n")
cat("Cell types:   ", paste(cell_types_keep, collapse=", "), "\n")
cat("Slot:         ", opt$slot, "\n")
cat("Inicio:       ", format(Sys.time()), "\n\n")

# ── 1. CARGAR OBJETOS ─────────────────────────────────────────────────────────
cat("[1/4] Cargando Seurat...\n")
obj <- readRDS(opt$seurat_rds)
cat(sprintf("      %s células | %s genes\n",
            format(ncol(obj), big.mark=","),
            format(nrow(obj), big.mark=",")))
#obj <- JoinLayers(obj)
#cat("      Layers unidos OK\n")


cat("[2/4] Cargando SC object...\n")
SC <- readRDS(opt$sc_rds)
cat(sprintf("      %s metacélulas | %s células\n",
            format(SC$N.SC, big.mark=","),
            format(length(SC$membership), big.mark=",")))

# ── 2. FILTRAR NAs ────────────────────────────────────────────────────────────
cat("[3/4] Filtrando células sin anotación...\n")
cells_keep <- !is.na(obj@meta.data[[opt$cell_type_col]]) &
              !is.na(obj@meta.data[[opt$phenotype_col]]) &
              obj@meta.data[[opt$cell_type_col]] != "" &
              obj@meta.data[[opt$phenotype_col]] != ""

obj_clean <- obj[, cells_keep]
cat(sprintf("      Células limpias: %s\n", format(sum(cells_keep), big.mark=",")))
#obj_clean <- JoinLayers(obj_clean)       # ← aquí
#cat("      Layers unidos OK\n") 
# ── 3. AGREGAR EXPRESION ──────────────────────────────────────────────────────
cat("[4/4] Agregando expresión con supercell_GE()...\n")

# Extraer solo layers data (no counts ni scale.data)
data_layers <- names(obj_clean[["RNA"]]@layers)
data_layers <- data_layers[grepl("^data\\.", data_layers)]

ge_matrix_full <- do.call(cbind, lapply(
  data_layers,
  function(l) obj_clean[["RNA"]]@layers[[l]]
))

ge_matrix <- supercell_GE(
  ge     = ge_matrix_full,
  groups = SC$membership
)

# Agregar nombres de genes y metacélulas
rownames(ge_matrix) <- rownames(ge_matrix_full)
colnames(ge_matrix) <- paste0("MC", 1:ncol(ge_matrix))

cat(sprintf("      Matriz agregada: %s genes × %s metacélulas\n",
            format(nrow(ge_matrix), big.mark=","),
            format(ncol(ge_matrix), big.mark=",")))
# ── 4. SPLIT Y GUARDAR TSVs ───────────────────────────────────────────────────
cat("\nGenerando matrices por cell type × fenotipo...\n")

cell_type_sc <- SC$SC.cell.annotation.
phenotype_sc <- SC$SC.cell.split.condition.

resumen  <- data.frame()
contador <- 0

for(ct in cell_types_keep){

  if(!ct %in% unique(cell_type_sc)){
    cat(sprintf("  [WARN] '%s' no encontrado en el objeto SC, saltando...\n", ct))
    next
  }

  for(ph in unique(phenotype_sc)){

    idx <- which(cell_type_sc == ct & phenotype_sc == ph)

    if(length(idx) == 0){
      cat(sprintf("  [SKIP] %s × %s — sin metacélulas\n", ct, ph))
      next
    }

    mat <- ge_matrix[, idx, drop = FALSE]

    ct_clean <- gsub(" ", "_", gsub("/", "_", ct))
    ph_clean <- gsub(" ", "_", ph)
    fname    <- file.path(opt$out_dir, paste0(ct_clean, "_", ph_clean, ".tsv"))

    write.table(mat,
                file      = fname,
                sep       = "\t",
                quote     = FALSE,
                row.names = TRUE,
                col.names = TRUE)

    contador <- contador + 1
    cat(sprintf("  [%02d] %s × %s → %s genes × %s metacélulas → %s\n",
                contador, ct, ph,
                format(nrow(mat), big.mark=","),
                ncol(mat),
                basename(fname)))

    resumen <- rbind(resumen, data.frame(
      cell_type   = ct,
      phenotype   = ph,
      n_metacells = ncol(mat),
      n_genes     = nrow(mat),
      file        = basename(fname)
    ))
  }
}

# ── 5. RESUMEN ────────────────────────────────────────────────────────────────
write.csv(resumen,
          file.path(opt$out_dir, "resumen_matrices.csv"),
          row.names = FALSE)

cat(sprintf("\n✓ %d matrices guardadas en: %s\n", contador, opt$out_dir))
cat(sprintf("✓ Resumen: %s\n", file.path(opt$out_dir, "resumen_matrices.csv")))
cat(sprintf("Fin: %s\n", format(Sys.time())))
print(resumen)
