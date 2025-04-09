#Script para llamar las %%MatrixMarket matrix coordinate integer general
#y la metadata para
#1. generar matrices de conteos por lote
#2. Generar matrices de conteos por condicion 
#3. Matrices por tipo celular

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

#Parallel implementation plan for future

future::plan(multisession, workers = 10)

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

demultiplex.f <- function(sample_name, seurat_list, demultiplex) {
  message("Processing batch: ", sample_name)
  
  #Pick a batch
  seurat_obj <- seurat_list[[sample_name]]
  
  #Filter metadata for this batch
  meta_batch <- demultiplex %>%
    filter(libraryBatch == sample_name)
  
  #Join with barcodes
  indiv_df <- tibble(cellBarcode = colnames(seurat_obj)) %>%
    left_join(meta_batch, by = "cellBarcode")%>%
    filter(!is.na(individualID))  #Filter NAs
  
  #Subset Seurat object
  seurat_obj <- subset(seurat_obj, cells = indiv_df$cellBarcode)
  
  #Assign individualID
  seurat_obj$individualID <- indiv_df$individualID
  
  return(seurat_list)
}

#Get data --- ---
cat("#Get data --- ---")

#Set directory
directory <- "/datos/rosmap/single_cell/matrix_exp_2-minimal/"

#Obtain matrix files .mtx.gz

matrix_files <- list.files(directory,
                           pattern = "\\.matrix\\.mtx\\.gz$", full.names = TRUE)

#Extract names from samples
batch_names <- str_replace(basename(matrix_files), "\\.matrix\\.mtx\\.gz$", "")

#List of Seurat objects
seurat_list <- future_map(batch_names, read_matrix)
names(seurat_list) <- batch_names

#Demultiplex --- ---

cat("#Demultiplex --- ---")

demultiplex <- vroom(file = "/datos/rosmap/single_cell/metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv")

seurat_demux_list <- future_map(
  sample_names,
  ~ demultiplex.f(.x, seurat_list, demultiplex)
)

#Assign names
names(seurat_demux_list) <- sample_names
