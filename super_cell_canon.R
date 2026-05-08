#The first ateempt to run super cell 


R.home()

library(Seurat)
library(SuperCell)


obj <- readRDS("/STORAGE/csbig/sc_ADers/celltypist/unassigned_only_minimal_2/annotation_with_celltypist/merged_harmony_integrated_annotated_plus_celltypist.rds")
dim(obj)

# Filtrar NAs y valores vacíos
cells_keep <- !is.na(obj@meta.data$cell_type) & 
  !is.na(obj@meta.data$is_AD) &
  obj@meta.data$cell_type != "" &
  obj@meta.data$is_AD != ""

obj_clean <- obj[, cells_keep]

# Extraer embedding de Harmony
harmony_embedding <- Embeddings(obj_clean, reduction = "harmony")

# Extraer anotaciones y dropear niveles vacíos de factores
cell_type <- droplevels(as.factor(obj_clean@meta.data$cell_type))
fenotipo  <- droplevels(as.factor(obj_clean@meta.data$is_AD))

# Verificar que no hay NAs ni vacíos
cat("NAs en cell_type:", sum(is.na(cell_type)), "\n")
cat("NAs en fenotipo: ", sum(is.na(fenotipo)), "\n")
cat("Niveles cell_type:", levels(cell_type), "\n")
cat("Niveles fenotipo: ", levels(fenotipo), "\n")


#Build metacells 

SC <- SCimplify_from_embedding(
  X                    = harmony_embedding,
  gamma                = 20, #graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)
  k.knn                = 5,# parameter to compute single-cell kNN network
  n.pc                 = 30, #number of principal components to use for construction of single-cell kNN network
  do.approx            = FALSE, # compute approximate kNN in case of a large dataset (>50'000) , maybe for the complete dataset 
  approx.N             = 50000, # number of cells to subsample for an approximate approach
  seed                 = 42,
  cell.annotation      = cell_type,
  cell.split.condition = fenotipo
)

#Results
cat("Células originales:", ncol(obj), "\n")
cat("Metacélulas:", SC$N.SC, "\n")
cat("Gamma usado:", SC$gamma, "\n")
# Ver distribución de metacélulas por tipo celular y fenotipo
table(SC$SC.cell.annotation., SC$SC.cell.split.condition.)

#Plots 
# Plot 1 - por tipo celular
png("supercell_by_celltype.png", width=1200, height=1000, res=150)
supercell_plot(SC$graph.supercells,
               group = SC$SC.cell.annotation.,
               main  = paste("Metacell network by cell type, gamma =", SC$gamma),
               seed  = 42)
dev.off()

# Plot 2 - por fenotipo
png("supercell_by_phenotype.png", width=1200, height=1000, res=150)
supercell_plot(SC$graph.supercells,
               group = SC$SC.cell.split.condition.,
               main  = paste("Metacell network by phenotype, gamma =", SC$gamma),
               seed  = 42)
dev.off()

cat("Plots guardados como PNG\n")
#Save 
# Guardar con fecha y gamma automáticos
fecha     <- format(Sys.Date(), "%Y-%m-%d")
gamma_val <- SC$gamma
out_dir   <- "/STORAGE/csbig/sc_ADers/supercell"

# Crear directorio si no existe
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Nombre del archivo automático
nombre <- file.path(out_dir, 
                    paste0("SC_metacells_minimal_gamma", gamma_val, "_", fecha, ".rds"))

saveRDS(SC, file = nombre)
cat("Objeto guardado en:", nombre, "\n")