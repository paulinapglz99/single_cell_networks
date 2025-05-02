#!/usr/bin/env Rscript

start.time <- Sys.time()

# Libraries
message("#Libraries --- ---")
pacman::p_load("tidyverse", "Seurat", "SingleCellExperiment", "furrr", "vroom", "optparse")

# Command-line options
option_list <- list(
  make_option(c("-d", "--directory"), type = "character", help = "Directory with MTX matrices", metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", help = "File with demultiplexing data", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", help = "Output file .rds", metavar = "character"),
  make_option(c("-t", "--test"), action = "store_true", default = FALSE,
              help = "Run test mode (only one batch)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$directory) | is.null(opt$metadata) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Mandatory arguments are missing: --directory, --metadata and/or --output", call. = FALSE)
}

# Parallel setup
message("#Parallel setup --- ---")
set.seed(10)
future::plan(multisession, workers = 10)

message("#Set functions --- ---")

# Demultiplexing function
demultiplex.f <- function(batch_name) {
  message("Reading and demultiplexing batch: ", batch_name)
  
  counts <- ReadMtx(
    mtx = file.path(opt$directory, paste0(batch_name, ".matrix.mtx.gz")),
    features = file.path(opt$directory, paste0(batch_name, ".features.tsv.gz")),
    cells = file.path(opt$directory, paste0(batch_name, ".barcodes.tsv.gz"))
  )
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = batch_name)
  seurat_obj$sample <- batch_name
  
  demux_subset <- demultiplex.df %>%
    filter(libraryBatch == batch_name)
  
  meta <- tibble(cell = colnames(seurat_obj)) %>%
    left_join(demux_subset, by = c("cell" = "cellBarcode")) %>%
    filter(!is.na(individualID)) %>%
    column_to_rownames("cell")
  
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
  seurat_ind <- SplitObject(seurat_obj, split.by = "individualID")
  names(seurat_ind) <- paste0(batch_name, "_", names(seurat_ind))
  
  message(Sys.time(), " - Processed individuals: ", length(seurat_ind))
  return(seurat_ind)
}

# Load demultiplex metadata
message("#Loading metadata --- ---")
demultiplex.df <- vroom(opt$metadata) %>% as_tibble()

# Get batch names
matrix_files <- list.files(opt$directory, pattern = "\\.matrix\\.mtx\\.gz$", full.names = TRUE)
batch_names <- str_replace(basename(matrix_files), "\\.matrix\\.mtx\\.gz$", "")
if (opt$test) {
  message("# Test mode active: running on one batch only --- ---")
  batch_names <- batch_names[1]
}

# Run demultiplexing
message("#Running demultiplexing in parallel --- ---")
list_seurat_ind <- future_map(
  batch_names,
  demultiplex.f,
  .options = furrr_options(
    seed = TRUE,
    globals = c("demultiplex.df","opt"),
    packages = c("dplyr", "Seurat", "tibble")
  ),
  .progress = TRUE
)
names(list_seurat_ind) <- batch_names

# Flatten and save
message("#Flattening and saving --- ---")
list_seurat_ind <- flatten(list_seurat_ind)
saveRDS(list_seurat_ind, file = opt$output,  compress = "gzip")

message("Running time:", Sys.time() - start.time, " min")
#END
