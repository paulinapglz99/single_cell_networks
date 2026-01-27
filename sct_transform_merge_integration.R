#!/usr/bin/env Rscript

#This script takes single cell qced data and transforms it with sctransform
#then merges data by individual
#and finally, integrates data and correct batch and other technical effects

#If any of this packages needed..
if (!requireNamespace("pacman", quietly = FALSE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = FALSE)) install.packages("optparse", repos = "https://cloud.r-project.org")
#install.packages("sctransform")
#BiocManager::install('glmGamPoi')

ok <- pacman::p_load("Seurat",
                     "tidyverse",
                     "sctransform",
                     "purrr",
                     "optparse",
                     "future",
                     "vroom")

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
  make_option(c("-o", "--out_dir"), type = "character", default = "results_sct_integration", help = "Output directory [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer", default = 4, help = "Parallel workers [default: %default]"),
  make_option(c("--seed"),type = "integer", default = 42, help = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

stopifnot(file.exists(opt$seurat_list),
  file.exists(opt$assay_metadata),
  file.exists(opt$clinical_metadata))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

#Parallelization

plan(multicore, workers = opt$workers)
options(future.globals.maxSize = 200 * 1024^3)

#Load data

#Just por ahora, comentar despues

opt$seurat_list <- "/STORAGE/csbig/sc_ADers/qc_out/matrices_demultiplexed_minimal-2025-11-05_15-09/seurat_list_filtered.rds"
opt$assay_metadata <- "/STORAGE/csbig/sc_ADers/metadata/ROSMAP_assay_scrnaSeq_metadata.csv"
opt$clinical_metadata <- "/STORAGE/csbig/sc_ADers/metadata/tables/clinical_stratified.csv"

message("Loading Seurat list...")
x <- readRDS(opt$seurat_list)
x <- x[1:5]

message("Loading metadata...")
assay_metadata.df    <- vroom::vroom(opt$assay_metadata) %>%  distinct()
clinical_metadata.df <- vroom::vroom(opt$clinical_metadata) %>% 
  filter(!is.na(is_AD)) %>% distinct()

#Add metadata to the Seurat object

assay_metadata.df <- assay_metadata.df %>%
  mutate(individualID = stringr::str_extract(
      specimenID,
      "R[0-9]+$"))

assay_metadata.df <- assay_metadata.df %>%
  filter(individualID %in% clinical_metadata.df$individualID)

#Metadata final por specimen
meta_by_specimen <- assay_metadata.df %>%
  select(
    specimenID,
    individualID,
    sequencingBatch,
    libraryBatch,
    platformLocation,
    RIN
  ) %>%
  left_join(
    clinical_metadata.df,
    by = "individualID"
  )

valid_specimen_ids <- intersect(names(x),
  meta_by_specimen$specimenID)

message("Keeping ", length(valid_specimen_ids),
        " / ", length(x),
        " Seurat objects")
#Filter only specimen ids with a valid 
x <- x[valid_specimen_ids]

#Add metadata


x <- setNames(
  lapply(names(x), function(spec_id) {
    
    seu <- x[[spec_id]]
    
    meta_row <- meta_by_specimen %>%
      filter(specimenID == spec_id)
    
    stopifnot(nrow(meta_row) == 1)
    
    seu$specimenID       <- meta_row$specimenID
    seu$individualID     <- meta_row$individualID
    seu$sequencingBatch  <- meta_row$sequencingBatch
    seu$libraryBatch     <- meta_row$libraryBatch
    seu$platformLocation <- meta_row$platformLocation
    seu$RIN              <- meta_row$RIN
    
    seu$is_AD   <- meta_row$is_AD
    seu$cogdx   <- meta_row$cogdx
    seu$ceradsc <- meta_row$ceradsc
    seu$braaksc <- meta_row$braaksc
    seu$pmi     <- meta_row$pmi
    
    seu
  }),
  names(x)
)

#FALTA: ver como agregar metadata clinica y lo del merge final???

#Now apply SCT transform

plan(multicore, workers = 16)
options(future.globals.maxSize = 200 * 1024^3)

x <- map(x, ~ SCTransform(
  .x,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  vst.flavor = "v2",
  return.only.var.genes = FALSE, #change if it's too lazy
  verbose = TRUE
))

#saveRDS(x, "x_after_SCT.rds")

#Check PCA

#merge w/0 integration or merging
seu.merged_pca <- merge(x[[1]], x[-1])

DefaultAssay(seu.merged_pca) <- "SCT"
VariableFeatures(seu.merged_pca) <- VariableFeatures(x[[1]])

seu.merged_pca <- RunPCA(
  seu.merged_pca,
  features = rownames(seu.merged_pca),
  npcs = 50,
  verbose = FALSE
)

#Vis UMAPs

pdf("pca_wo_integration_platformLocation.pdf")
DimPlot(seu.merged_pca, group.by = "platformLocation")
dev.off()

pdf("pca_wo_integration_libraryBatch.pdf")
DimPlot(seu.merged_pca, group.by = "libraryBatch")
dev.off()

#Next steps are

#Select Integration Features

features <- SelectIntegrationFeatures(
  object.list = x,
  nfeatures = 3000)

# for (obj.i in merged_by_individual) {
#   sct_genes <- rownames(obj.i[["SCT"]]@scale.data)
#   features <- intersect(features, sct_genes)
# }

#PrepSCTIntegration

x <- PrepSCTIntegration(
  object.list = x,
  anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = x,
  normalization.method = "SCT",
  anchor.features = features)

seu.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT")

seu.integrated <- RunPCA(seu.integrated, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
seu.integrated <- FindNeighbors(seu.integrated, reduction = "pca", dims = 1:30)
seu.integrated <- FindClusters(seu.integrated, resolution = 0.3)

#PCA

pdf("pca_wo_integration_merged_platformLocation.pdf")
DimPlot(seu.integrated, reduction = "umap",
        #split.by = "stim",
        group.by = "libraryBatch")
dev.off()

#Finally 

#Merge Seurat objects by individual

#For each object, we take the unique individualID
individual_per_object <- sapply(seu.integrated, function(seu) {
  unique(seu@meta.data$individualID)
})

#Sanity check: each object must have only one individual :$
stopifnot(all(lengths(individual_per_object) == 1))

#Group objects by individual, this creates a list where each element is
#a list of Seurat objects from the same individual
objects_ind <- split(x, individual_per_object)

#Merge libraries belonging to the same individual

merged_by_individual <- lapply(
  objects_ind,
  function(seu_list) {
    
    #If there is only one bookshop, return as is.
    if (length(seu_list) == 1) {
      return(seu_list[[1]])
    }
    
    #Initialize for loop
    merged_seu <- seu_list[[1]]
    
    #Merge everything
    for (i in 2:length(seu_list)) {
      merged_seu <- merge(
        merged_seu,
        seu_list[[i]],
        merge.data = TRUE
      )}
    
    merged_seu
  })

#Check merged_by_individual. It's a list of Seurat objs, one per individual
length(merged_by_individual)

#Vis PCA again

seu.merged_pca_postmerge <- merge(seu.integrated[[1]], seu.integrated[-1])

DefaultAssay(seu.merged_pca_postmerge) <- "SCT"
VariableFeatures(seu.merged_pca_postmerge) <- VariableFeatures(x[[1]])

seu.merged_pca_postmerge <- RunPCA(
  seu.merged_pca_postmerge,
  features = rownames(seu.merged_pca_postmerge),
  npcs = 50,
  verbose = FALSE
)

#Vis UMAPs

pdf("pca_wo_integration_merged_platformLocation.pdf")
DimPlot(seu.merged_pca, group.by = "platformLocation")
dev.off()

pdf("pca_wo_integration_merged_libraryBatch.pdf")
DimPlot(seu.merged_pca, group.by = "libraryBatch")
dev.off()

#Save output
output_file <- file.path(outdir, "seurat_list_sctransform_normalized.rds")
saveRDS(seurat_list_norm, file = output_file)
message("#DONE. Saved in: ", output_file)