#Script to call for cell matrix
#and metadata for
#Demultiplexing count matrices per individual

#Libraries --- ---
cat("#Libraries --- ---")

#install.packages('Seurat')

pacman::p_load("tidyverse", 
               "Seurat", 
               "ggplot2", 
               "SingleCellExperiment", 
               "furrr", 
               "vroom")

#Load functions --- ---
cat("#Load functions --- ---")

#Parallel implementation plan for future
set.seed(10)
future::plan(multisession, workers = 5)

#Function to read matrix

read_matrix.f <- function(sample) {
  message("Reading sample: ", sample)
  
  counts <- ReadMtx(
    mtx = file.path(directory, paste0(sample, ".matrix.mtx.gz")),
    features = file.path(directory, paste0(sample, ".features.tsv.gz")),
    cells = file.path(directory, paste0(sample, ".barcodes.tsv.gz"))
  )
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  seurat_obj$sample <- sample
  return(seurat_obj)
}


#Function to demultiplex matrix

demultiplex.f <- function(batch_name) {
  #Filter only corresponding batch
  demux_subset <- demultiplex.df %>% filter(libraryBatch == batch_name)
  #Cell barcodes of the batch
  barcodes <- names(seurat_list[[batch_name]]$sample)
  #Match demultiplex data with barcodes
  meta <- data.frame(cell = barcodes) %>%
    left_join(demux_subset, by = c("cell" = "cellBarcode")) %>%
    filter(!is.na(individualID)) %>% 
    column_to_rownames("cell")
  #Assign metadata to the Seurat object
  seurat_obj <- AddMetaData(seurat_list[[batch_name]], metadata = meta)
  seurat_ind <- SplitObject(seurat_obj, split.by = "individualID")
  #Rename element in format batch_individualID
  names(seurat_ind) <- paste0(batch_name, "_", names(seurat_ind))
  
  message(Sys.time(), " - Processed individuals: ", length(seurat_ind))
  return(seurat_ind)
  
}

#Get data --- ---
cat("#Get data --- ---")

#Set directory
directory <- "/datos/rosmap/single_cell/matrix_exp_2-minimal/"

#Obtain matrix files .mtx.gz

matrix_files <- list.files(directory,
                           pattern = "\\.matrix\\.mtx\\.gz$", full.names = TRUE)

#Extract names from batch samples

batch_names <- names(seurat_list)

#List of Seurat objects
seurat_list <- future_map(batch_names, read_matrix.f)
names(seurat_list) <- batch_names

#Demultiplex --- ---

#Try this function 
# 
# seurat_obj <- seurat_list[1]
# batch_name <- batch_names[1]
# 
#test <- demultiplex.f(batch_name = batch_names[1])

#Demultiplex --- ---

cat("#Demultiplex --- ---")

demultiplex.df <- vroom(file = "/datos/rosmap/single_cell/metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv")

demux_list <- split(demultiplex.df, demultiplex.df$libraryBatch)

#Parallel version using future_map

resultados <- future_map(
  batch_names, 
  demultiplex.f,
  .options = furrr_options(
    #scheduling = 2,
    seed = TRUE,           # Para reproducibilidad
    globals = c("seurat_list", "demultiplex.df"),  # Exportar variables a los workers
    packages = c("dplyr", "Seurat", "tibble")  # Paquetes necesarios en workers
  ),
  .progress = TRUE        # Barra de progreso
)

