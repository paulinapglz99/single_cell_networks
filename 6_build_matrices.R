#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Matrix)
})

# ── PARÁMETROS ──────────────────────────────────────────────────────────────
seurat_rds  <- "/STORAGE/csbig/sc_ADers/celltypist/celltypist/unassigned_only_final/annotation_with_celltypist/merged_harmony_integrated_annotated_plus_celltypist.rds"
sc_rds      <- "/STORAGE/csbig/sc_ADers/supercell/SC_gamma20_2026-05-08.rds"
out_dir     <- "/STORAGE/csbig/sc_ADers/supercell/matrices_aracne"
cell_type_col <- "cell_type"
phenotype_col <- "is_AD"

cell_types_keep <- c(
  "Excitatory Neurons",
  "Inhibitory Neurons",
  "Astrocyte",
  "Oligodendrocytes",
  "Endothelial",
  "OPCs",
  "Microglia"
)

# ── SETUP ────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Salida:", out_dir, "\n")

# ── 1. CARGAR OBJETOS ────────────────────────────────────────────────────────
cat("[1/4] Cargando Seurat...\n")
obj <- readRDS(seurat_rds)
cat(sprintf("      %s células | %s genes\n",
            format(ncol(obj), big.mark=","),
            format(nrow(obj), big.mark=",")))

cat("[2/4] Cargando SC gamma20...\n")
SC <- readRDS(sc_rds)
cat(sprintf("      %s metacélulas | %s células\n",
            format(SC$N.SC, big.mark=","),
            format(length(SC$membership), big.mark=",")))

# ── 2. FILTRAR NAs (igual que en supercell script) ───────────────────────────
cat("[3/4] Filtrando células sin anotación...\n")
cells_keep <- !is.na(obj@meta.data[[cell_type_col]]) &
              !is.na(obj@meta.data[[phenotype_col]]) &
              obj@meta.data[[cell_type_col]] != "" &
              obj@meta.data[[phenotype_col]] != ""

obj_clean <- obj[, cells_keep]
cat(sprintf("      Células limpias: %s\n", format(sum(cells_keep), big.mark=",")))

# ── 3. AGREGAR EXPRESIÓN ─────────────────────────────────────────────────────
cat("[4/4] Agregando expresión con supercell_GE()...\n")
ge_matrix <- supercell_GE(
  ge     = GetAssayData(obj_clean, slot = "data"),
  groups = SC$membership
)
cat(sprintf("      Matriz agregada: %s genes × %s metacélulas\n",
            format(nrow(ge_matrix), big.mark=","),
            format(ncol(ge_matrix), big.mark=",")))

# ── 4. SPLIT POR CELL TYPE × FENOTIPO Y GUARDAR TSV ─────────────────────────
cat("\nGenerando matrices por cell type × fenotipo...\n")

cell_type_sc <- SC$SC.cell.annotation.
phenotype_sc <- SC$SC.cell.split.condition.

resumen <- data.frame()
contador <- 0

for(ct in cell_types_keep){
  for(ph in unique(phenotype_sc)){
    
    # índices de metacélulas que pertenecen a este ct × fenotipo
    idx <- which(cell_type_sc == ct & phenotype_sc == ph)
    
    if(length(idx) == 0){
      cat(sprintf("  [SKIP] %s × %s — sin metacélulas\n", ct, ph))
      next
    }
    
    # submatriz
    mat <- ge_matrix[, idx, drop = FALSE]
    
    # nombre de archivo limpio
    ct_clean <- gsub(" ", "_", gsub("/", "_", ct))
    ph_clean <- gsub(" ", "_", ph)
    fname    <- file.path(out_dir, paste0(ct_clean, "_", ph_clean, ".tsv"))
    
    # guardar con genes como rownames
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

# ── 5. RESUMEN ───────────────────────────────────────────────────────────────
write.csv(resumen,
          file.path(out_dir, "resumen_matrices.csv"),
          row.names = FALSE)

cat(sprintf("\n✓ %d matrices guardadas en: %s\n", contador, out_dir))
cat("✓ Resumen guardado en: resumen_matrices.csv\n")
print(resumen)
