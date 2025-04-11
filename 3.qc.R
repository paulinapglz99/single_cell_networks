#!/usr/bin/env Rscript

#Script fot quality control of single cell data

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

#Load functions --- ---
message("#Load functions --- ---")

#Parallel implementation plan for future
set.seed(10)
future::plan(multisession, workers = 5)

#Get data --- ---
message("#Get data --- ---")

seurat_ind <- readRDS(file = "/datos/rosmap/single_cell/matrices_demultiplexed.rds")

#Mitocondrial % per sample
message("#Assign mitocondrial % per sample --- ---")

seurat_ind <- lapply(seurat_ind, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")  # si son genes humanos
  return(obj)
})

#Thresholds

min_genes <- 200
max_genes <- 6000
max_mt <- 10



