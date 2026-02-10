#!/usr/bin/env Rscript

#This script takes single cell qced data and transforms it, 
#then merges data by individual
#and finally, integrates data and correct batch and other technical effects

start_time <- Sys.time()
#If any of this packages needed..
if (!requireNamespace("pacman", quietly = FALSE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = FALSE)) install.packages("optparse", repos = "https://cloud.r-project.org")
#install.packages("sctransform")
#BiocManager::install('glmGamPoi')

ok <- pacman::p_load("Seurat",
                     "tidyverse",
                     #"sctransform",
                     "purrr",
                     "optparse",
                     "future",
                     "vroom", 
                     "harmony")

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages loaded correctly.: ",
       paste(names(ok)[!ok], collapse = ", "))
}

#Define option list for inputs

option_list <- list(
  make_option(c("-s", "--seurat_list"), type = "character", help = "Path to QCed Seurat list (.rds)"),
  make_option(c("-a", "--assay_metadata"),type = "character", help = "Assay metadata CSV"),
  make_option(c("-c", "--clinical_metadata"), type = "character", help = "Clinical metadata CSV"),
  make_option(c("-o", "--out_dir"), type = "character", default = "results_merge_integration", help = "Output directory [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer", default = 4, help = "Parallel workers [default: %default]"),
  make_option(c("--seed"),type = "integer", default = 42, help = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

stopifnot(file.exists(opt$seurat_list),
          file.exists(opt$assay_metadata),
          file.exists(opt$clinical_metadata))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(opt$out_dir)

#Parallelization

plan(multicore, workers = opt$workers)
options(future.globals.maxSize = 200 * 1024^3)

#Load data

#Just por ahora, comentar despues
# 
# opt$seurat_list <- "/STORAGE/csbig/matrices_demultiplexed_minimal_QC-2026-01-28_19-04/seurat_list_filtered.rds"
# opt$assay_metadata <- "/STORAGE/csbig/sc_ADers/metadata/ROSMAP_assay_scrnaSeq_metadata.csv"
# opt$clinical_metadata <- "/STORAGE/csbig/sc_ADers/metadata/tables/clinical_stratified.csv"
# opt$out_dir <- "~/AD_single_cell/merge_integration_minimal_data_frito"

message("Loading Seurat list...")
x <- readRDS(opt$seurat_list)
#x <- x[1:5]

message("Loading metadata...")
assay_metadata.df    <- vroom::vroom(opt$assay_metadata) %>%  distinct()
clinical_metadata.df <- vroom::vroom(opt$clinical_metadata) %>% 
  filter(!is.na(is_AD)) %>% distinct()
dim(assay_metadata.df)
dim(clinical_metadata.df)
#Add metadata to the Seurat object

assay_metadata.df <- assay_metadata.df %>%
  mutate(individualID = stringr::str_extract(
    specimenID,
    "R[0-9]+$"))

assay_metadata.df <- assay_metadata.df %>%
  filter(individualID %in% clinical_metadata.df$individualID)
dim(assay_metadata.df)
#Metadata final por specimen
meta_by_specimen <- assay_metadata.df %>%
  select(
    specimenID,
    individualID,
    sequencingBatch,
    libraryBatch,
    platformLocation,
    RIN) %>%
  left_join(
    clinical_metadata.df,
    by = "individualID") %>%
  mutate(specimenID_Location = paste(specimenID, platformLocation, sep = "_"))

valid_specimen_ids <- intersect(names(x),
                                meta_by_specimen$specimenID)

message("Keeping ", length(valid_specimen_ids),
        " / ", length(x),
        " Seurat objects")
#Filter only specimen ids with a valid 
x <- x[valid_specimen_ids]

name_map <- meta_by_specimen %>%
  distinct(specimenID, platformLocation, specimenID_Location)

new_names <- sapply(names(x), function(old_name) {
  rows <- name_map %>% filter(specimenID == old_name)
  
  if (nrow(rows) != 1) {
    stop("Ambiguous specimenID: ", old_name)
  }
  
  rows$specimenID_Location
})

names(x) <- new_names

message("#Add metadata \n")

x <- setNames(
  lapply(names(x), function(id) {
    
    seu <- x[[id]]
    
    meta_row <- meta_by_specimen %>%
      filter(specimenID_Location == id)
    
    stopifnot(nrow(meta_row) == 1)
    
    seu$specimenID            <- meta_row$specimenID
    seu$specimenID_Location   <- meta_row$specimenID_Location
    seu$individualID          <- meta_row$individualID
    seu$sequencingBatch       <- meta_row$sequencingBatch
    seu$libraryBatch          <- meta_row$libraryBatch
    seu$platformLocation      <- meta_row$platformLocation
    seu$RIN                   <- meta_row$RIN
    
    seu$is_AD   <- meta_row$is_AD
    seu$cogdx   <- meta_row$cogdx
    seu$ceradsc <- meta_row$ceradsc
    seu$braaksc <- meta_row$braaksc
    seu$pmi     <- meta_row$pmi
    
    seu
  }),
  names(x)
)

# #Merge Seurat objects by individual
# 
# #For each object, we take the unique individualID
# individual_per_object <- sapply(x, function(seu) {
#   unique(seu@meta.data$individualID)
# })
# 
# #Sanity check: each object must have only one individual :$
# stopifnot(all(lengths(individual_per_object) == 1))
# 
# #Group objects by individual, this creates a list where each element is
# #a list of Seurat objects from the same individual
# objects_ind <- split(x, individual_per_object)
# 
# #Merge libraries belonging to the same individual
# 
# merged_by_individual <- lapply(
#   objects_ind,
#   function(seu_list) {
#     
#     #If there is only one bookshop, return as is.
#     if (length(seu_list) == 1) {
#       return(seu_list[[1]])
#     }
#     
#     #Initialize for loop
#     merged_seu <- seu_list[[1]]
#     
#     #Merge everything
#     for (i in 2:length(seu_list)) {
#       merged_seu <- merge(
#         merged_seu,
#         seu_list[[i]],
#         merge.data = TRUE
#       )}
#     
#     merged_seu
#   })
# 
# #Check merged_by_individual. It's a list of Seurat objs, one per individual
# length(merged_by_individual)

#Merge everything lol
merged <- merge(x[[1]], x[-1])

#Normalization

DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged, verbose = FALSE)

#Find Variable Features

merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

#Scale data

merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)

#Run PCA

merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = 50, verbose = FALSE)

#Check elbowplot
pdf("ElbowPlot_merged.pdf")
ElbowPlot(merged, ndims = 50)
dev.off()

merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)

#Vis UMAPs

pdf(file.path(opt$out_dir,"umap_merged_individual_sequencingBatch_preintegration.pdf"))
DimPlot(merged,
        reduction = "umap",
        group.by = "sequencingBatch")
dev.off()

pdf(file.path(opt$out_dir,"umap_merged_individual_platformLocation_preintegration.pdf"))
DimPlot(merged,
        reduction = "umap",
        group.by = "platformLocation")
dev.off()

pdf(file.path(opt$out_dir,"umap_merged_individual_is_AD_preintegration.pdf"))
DimPlot(merged,
        reduction = "umap",
        group.by = "is_AD")
dev.off()

message("#Run harmony and check for batches \n")

merged$sequencingBatch  <- factor(merged$sequencingBatch)
merged$platformLocation <- factor(merged$platformLocation)

batch_vars <- c("sequencingBatch"#,
                #              "platformLocation"
)

merged <- RunHarmony(
  merged,
  group.by.vars = batch_vars,
  #reduction = "pca",
  assay.use = "RNA",
  verbose = TRUE)

message("#Rum UMAP in new harmonized data \n")

merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:30)

message("#Find FindNeighbors \n")

merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims = 1:30)

merged <- FindClusters(
  merged,
  resolution = 0.3)

#Vis UMAP
pdf(file.path(opt$out_dir,"umap_harmony_sequencingBatch.pdf"))
DimPlot(merged, reduction = "umap", group.by = "sequencingBatch")
dev.off()

pdf(file.path(opt$out_dir,"umap_harmony_platformLocation.pdf"))
DimPlot(merged, reduction = "umap", group.by = "platformLocation")
dev.off()

pdf(file.path(opt$out_dir,"umap_harmony_is_AD.pdf"))
DimPlot(merged, reduction = "umap", split.by = "is_AD")
dev.off()

#Finally 

message("#Save output \n")
saveRDS(merged, file = file.path(opt$out_dir, "merged_by_individual_harmony.rds"))

message("#Save metadata\n")
write.csv(
  merged@meta.data,
  file = file.path(opt$out_dir, "merged_by_individual_harmony_cell_metadata.csv"),
  row.names = TRUE)

message("#DONE. Saved in: ", opt$out_dir)
end_time <- Sys.time()
time_taken <- end_time - start_time
message(print(time_taken))

message("#THE END?\n")
