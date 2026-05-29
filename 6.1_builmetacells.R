# Build metacells based on 5.supercell_graph
# Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Matrix)
  library(optparse)
})

# Arguments
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
              help    = "Comma-separated cell types [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validation
stopifnot(file.exists(opt$seurat_rds))
stopifnot(file.exists(opt$sc_rds))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Load objects
obj <- readRDS(opt$seurat_rds)
SC  <- readRDS(opt$sc_rds)

# 1. Filtrar NAs primero
cells_keep <- !is.na(obj@meta.data[[opt$cell_type_col]]) &
  !is.na(obj@meta.data[[opt$phenotype_col]]) &
  obj@meta.data[[opt$cell_type_col]] != "" &
  obj@meta.data[[opt$phenotype_col]] != ""

obj_clean <- obj[, cells_keep]

# 1.2 JoinLayers DESPUÉS del filtro
#obj_clean <- JoinLayers(obj_clean)

# 2. Extract expression matrix -Esta linea es la buena 
#ge_matrix_cells <- GetAssayData(obj_clean, assay = "RNA", layer = "data")

#probemos extarer matrix sIN S5

data_layers <- grep("^data\\.", names(obj_clean[["RNA"]]@layers), value = TRUE)

ge_matrix_cells <- do.call(cbind, lapply(
  data_layers,
  function(l) obj_clean[["RNA"]]@layers[[l]]
))

rownames(ge_matrix_cells) <- rownames(obj_clean[["RNA"]])



# 3. Average with supercell_GE — principal function
ge_matrix <- supercell_GE(
  ge     = ge_matrix_cells,
  groups = SC$membership
)

rownames(ge_matrix) <- rownames(ge_matrix_cells)
colnames(ge_matrix) <- paste0("MC", seq_len(ncol(ge_matrix)))

# 4. Split by cell type × phenotype and save TSVs
cell_types_keep <- trimws(strsplit(opt$cell_types, ",")[[1]])
cell_type_sc    <- SC$SC.cell.annotation.
phenotype_sc    <- SC$SC.cell.split.condition.

resumen  <- data.frame()
contador <- 0

for (ct in cell_types_keep) {
  
  if (!ct %in% unique(cell_type_sc)) {
    cat(sprintf("  [WARN] '%s' no encontrado en SC, saltando...\n", ct))
    next
  }
  
  for (ph in unique(phenotype_sc)) {
    
    idx <- which(cell_type_sc == ct & phenotype_sc == ph)
    
    if (length(idx) == 0) {
      cat(sprintf("  [SKIP] %s x %s — sin metacélulas\n", ct, ph))
      next
    }
    
    mat <- ge_matrix[, idx, drop = FALSE]
    
    ct_clean <- gsub("[ /]", "_", ct)
    ph_clean <- gsub(" ", "_", ph)
    fname    <- file.path(opt$out_dir, paste0(ct_clean, "_", ph_clean, ".tsv"))
    
    write.table(mat,
                file      = fname,
                sep       = "\t",
                quote     = FALSE,
                row.names = TRUE,
                col.names = TRUE)
    
    contador <- contador + 1
    cat(sprintf("  [%02d] %s x %s → %s genes x %s MCs → %s\n",
                contador, ct, ph,
                format(nrow(mat), big.mark = ","),
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

# 5. Summary
write.csv(resumen,
          file.path(opt$out_dir, "resumen_matrices.csv"),
          row.names = FALSE)

cat(sprintf("\n✓ %d matrices guardadas en: %s\n", contador, opt$out_dir))
cat(sprintf("Fin: %s\n", format(Sys.time())))
print(resumen)

