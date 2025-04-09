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

future::plan(multisession, workers = 5)

#Function to read matrix

read_matrix <- function(sample) {
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

demultiplex.f <- function(seurat_obj, batch_name, demux_map) {
  #Original barcodes
  barcodes <- colnames(seurat_obj)
  
  # Filtrar solo el batch correspondiente
  demux_subset <- demux_map %>% filter(libraryBatch == batch_name)
  
  # Empatar barcodes con demux info
  meta <- tibble(cell = barcodes) %>%
    left_join(demux_subset, by = c("cell" = "cellBarcode")) %>%
    mutate(
      individualID = ifelse(is.na(individualID), "Unknown", individualID),
      unique_cell_name = paste(batch_name, individualID, cell, sep = "_")
    )
  
  # Agregar metadata y renombrar celdas
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta %>% column_to_rownames("cell"))
  colnames(seurat_obj) <- meta$unique_cell_name
  
  return(seurat_obj)
}

#Get data --- ---
cat("#Get data --- ---")

#Set directory
directory <- "/datos/rosmap/single_cell/matrix_exp_2-minimal/"

#Obtain matrix files .mtx.gz

matrix_files <- list.files(directory,
                           pattern = "\\.matrix\\.mtx\\.gz$", full.names = TRUE)

#Extract names from batch samples
batch_names <- str_replace(basename(matrix_files), "\\.matrix\\.mtx\\.gz$", "")

#List of Seurat objects
seurat_list <- future_map(batch_names, read_matrix)
names(seurat_list) <- batch_names

#Demultiplex --- ---

cat("#Demultiplex --- ---")

demultiplex.df <- vroom(file = "/datos/rosmap/single_cell/metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv")

seurat_test <- demultiplex.f(
  seurat_obj = seurat_list[["190403-B4-A"]],
  batch_name = "190403-B4-A",
  demux_map = demultiplex.df
)







