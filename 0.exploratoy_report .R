#Data description

#Script to describe single cell data

start.time <- Sys.time()

# Set seed for reproducibility
set.seed(123)

#Libraries --- ---
message("#Libraries --- ---")

#BiocManager::install("scDblFinder")

pacman::p_load("tidyverse", 
               "Seurat", 
               "ggplot2", 
               "SingleCellExperiment", 
               "furrr", 
               "vroom", 
               "optparse", 
               "cowplot",
               "scDblFinder")

seurat_ind <- readRDS(file = "/datos/rosmap/single_cell/matrices_demultiplexed.rds")

#Functions ---- ---

get_seurat_summary <- function(sample_names) {
  summary_list <- lapply(sample_names, function(sample) {
    
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
  
  #Summary statistics
  total_cells <- sum(summary$num_cells)
  avg_cells <- mean(summary$num_cells)
  sd_cells <- sd(summary$num_cells)
  genes <- unique(summary$num_genes) # Asumiendo que todos tienen el mismo número de genes
  
  #Print message
  message("\nStatistic summary:\n")
  message("------------------------------------\n")
  message("Totall cells in all samples:", total_cells, "\n")
  message("Average cells per individual:", round(avg_cells, 1), "\n")
  message("Desviación estándar de células por individuo:", round(sd_cells, 1), "\n")
  message("Número de genes:", genes[1], "\n") # Tomamos el primero ya que deberían ser iguales
  message("------------------------------------\n\n")
  
  #Give df
  return(summary)
  
}

#How many cells and genes do we have? --- ---

sample_names <- names(seurat_ind)

summary <- get_seurat_summary(sample_names = sample_names)

#Calculate the percentage of mitochondrial genes --- ---

seurat_ind[sample_names] <- lapply(seurat_ind[sample_names], function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
})

seurat_ind[[1]]@meta.data %>% head()

vlPlot <- VlnPlot(seurat_ind[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
vlPlot

#

sample1.qc <- FetchData(seurat_ind[[1]],
                        vars=c("nFeature_RNA","nCount_RNA","percent.mt"))
nCountRNA <- sample1.qc %>%
  ggplot() +
  geom_histogram(aes(x=nCount_RNA), bins=100) +
  theme_minimal()
nFeatureRNA <- sample1.qc %>%
  ggplot() +
  geom_histogram(aes(x=nFeature_RNA), bins=100) +  
  theme_cowplot()

mt.count <- sample1.qc %>%
  ggplot() +
  geom_histogram(aes(x=percent.mt), bins=100) + 
  theme_cowplot()

nCountRNA + nFeatureRNA+mt.count

#Set an upper threshold of 6000 for nFeature_RNA and 35000 for *nCount_RNA

scPlot <- sample1.qc %>%
  mutate(keep = if_else(nCount_RNA > 35000 & nFeature_RNA > 6500, "remove", "keep")) %>%
  ggplot() +
  geom_point(aes(nCount_RNA, nFeature_RNA, colour=keep), alpha=.50) +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot()
scPlot

#Now the same for the lower treshold using  300 for nFeature_RNA** and 600 for nCount_RNA

scPlot <- sample1.qc %>%
  mutate(keep = if_else(nCount_RNA > 300 & nFeature_RNA > 600, "keep", "remove")) %>%
  ggplot() +
  geom_point(aes(nCount_RNA, nFeature_RNA, colour=keep), alpha=.50) +
  scale_x_log10() +
  scale_y_log10()
scPlot

#Identify and Remove Doublets
sce <- as.SingleCellExperiment(seurat_ind[[1]])
# Display the converted object
sce

## We need to set.seed() because the scDblFinder command
## uses a probabilist strategy to identify doublets. This
## means that, everytime we run the command, it will produce
## results that are slightly different. The set.seed()
## comman will guarantee the same results everytime.

db.results <- scDblFinder(sce, returnType = 'table') %>%
  as.data.frame() %>%
  filter(type == 'real')
head(db.results)


db.results %>%
  dplyr::count(db.results)

