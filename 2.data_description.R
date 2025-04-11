#Data description

#Script to describe single cell data

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

seurat_ind <- readRDS(file = "/datos/rosmap/single_cell/matrices_demultiplexed.rds")

#Functions ---- ---

get_seurat_summary <- function(seurat_list) {
  summary_df <- lapply(names(seurat_list), function(name) {
    obj <- seurat_list[[name]]  # Acceder al objeto correcto usando el nombre
    data.frame(
      sample = name,
      n_genes = nrow(obj@assays$RNA$counts),  # Genes = filas de la matriz
      n_cells = ncol(obj@assays$RNA$counts)   # Células = columnas de la matriz
    )
  }) %>% bind_rows()
  
  # Totales
  total_samples <- nrow(summary_df)
  total_genes <- sum(summary_df$n_genes)      # Suma de genes POR MUESTRA (no únicos)
  total_cells <- sum(summary_df$n_cells)
  
  # Mensaje resumen
  message("Summary of the dataset:")
  message("Individual samples: ", total_samples)
  message("Genes (suma por muestra): ", total_genes)
  message("Células (total): ", total_cells)
  
  return(summary_df)
}

#How many cells and genes do we have? --- ---

total_cells <- sum(sapply(seurat_ind, 
                          function(obj) ncol(obj[["RNA"]]$counts)))
total_genes <- length(unique(unlist(lapply(seurat_ind, function(obj) rownames(obj[["RNA"]]$counts)))))

message("The dataset contains ", total_genes, " genes", length(seurat_ind), " individuals and ",  total_cells, " cells.\n")

summary_seurat_list <- function(seurat_list) {
  stats <- lapply(seurat_list, function(obj) {
    c(Genes = nrow(obj@assays$RNA$counts), 
      Cells = ncol(obj@assays$RNA$counts))
  })
  
  total_genes <- sum(sapply(stats, `[`, "Genes"))
  total_cells <- sum(sapply(stats, `[`, "Cells"))
  
  message("Resumen del dataset:\n",
          "Muestras analizadas: ", length(seurat_list), "\n",
          "Genes totales (suma por muestra): ", total_genes, "\n",
          "Células totales: ", total_cells)
}



