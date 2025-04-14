#!/usr/bin/env Rscript

#Script to call for cell matrices and metadata for
#Demultiplexing count matrices per individual

start.time <- Sys.time()

#Libraries --- ---
message("#Libraries --- ---")

#install.packages('Seurat')

pacman::p_load("tidyverse", 
               "Seurat", 
               "ggplot2", 
               "SingleCellExperiment", 
               "furrr", 
               "vroom", 
               "optparse")

#Options from console --- ---

option_list <- list(
  make_option(c("-d", "--directory"), type = "character", help = "Directory with MTX matrices", metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", help = "File with demultiplexing data", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", help = "Output file .rds", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$directory) | is.null(opt$metadata) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Mandatory arguments are missing: --directory, --metadata and/or --output", call. = FALSE)
}

#Load functions --- ---
message("#Load functions --- ---")

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
message("#Get data --- ---")

#Set directory
#directory <- "/datos/rosmap/single_cell/matrix_exp_2-minimal/"
directory <- opt$directory

#Obtain matrix files .mtx.gz

matrix_files <- list.files(directory,
                           pattern = "\\.matrix\\.mtx\\.gz$", full.names = TRUE)

#Obtain batch names
batch_names <- str_replace(basename(matrix_files), "\\.matrix\\.mtx\\.gz$", "")

#List of Seurat objects
seurat_list <- future_map(batch_names, read_matrix.f)
names(seurat_list) <- batch_names

#Demultiplex --- ---

message("#Demultiplex --- ---")

#Try this function 
# 
# seurat_obj <- seurat_list[1]
# batch_name <- batch_names[1]
# 
#test <- demultiplex.f(batch_name = batch_names[1])

demultiplex.df <- vroom(opt$metadata)
#demultiplex.df <- vroom(file = "/datos/rosmap/single_cell/metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv")

#Parallel version using future_map

list_seurat_ind <- future_map(
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
names(list_seurat_ind) <- batch_names

#Flatten
message("#Flatten demultiplex --- ---")

list_seurat_ind <- flatten(list_seurat_ind)

#Save list

#saveRDS(list_seurat_ind, file = "/datos/rosmap/single_cell/matrices_demultiplexed.rds")
saveRDS(list_seurat_ind, file = opt$output)

message("Running time:", Sys.time() - start.time, " s")
#END