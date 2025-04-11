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

get_seurat_summary <- function(sample_names) {
  summary_list <- lapply(sample_names, function(sample) {
      
    # Si no está vacío, calcular normalmente
    num_cells <- length(seurat_ind[[sample]]@assays$RNA@cells)
    num_genes <- length(seurat_ind[[sample]]@assays$RNA@features)
    
    data.frame(
      sample = sample,
      num_cells = num_cells,
      num_genes = num_genes,
      stringsAsFactors = FALSE
    )
  })
  
  summary <- do.call(rbind, summary_list)
  
  # Calcular estadísticas summary
  total_cells <- sum(summary$num_cells)
  avg_cells <- mean(summary$num_cells)
  sd_cells <- sd(summary$num_cells)
  genes <- unique(summary$num_genes) # Asumiendo que todos tienen el mismo número de genes
  
  # Mostrar el mensaje con las estadísticas
  cat("\nResumen estadístico:\n")
  cat("------------------------------------\n")
  cat("Número total de células:", total_cells, "\n")
  cat("Promedio de células por individuo:", round(avg_cells, 1), "\n")
  cat("Desviación estándar de células por individuo:", round(sd_cells, 1), "\n")
  cat("Número de genes:", genes[1], "\n") # Tomamos el primero ya que deberían ser iguales
  cat("------------------------------------\n\n")
  
  # Devolver el dataframe como antes
  return(summary)
  
}


#How many cells and genes do we have? --- ---

sample_names <- names(seurat_ind)

summary <- get_seurat_summary(sample_names = sample_names)

