#!/bin/bash

# Test con muestra MINIMAL para aislar problema de memoria vs código
Rscript ~/single_cell_networks/3.merge_integration.R \
  -s /STORAGE/csbig/sc_ADers/matrices_demultiplexed_minimal_QC-2026-01-28_19-04/seurat_list_filtered.rds \
  -a /STORAGE/csbig/sc_ADers/metadata/ROSMAP_assay_scrnaSeq_metadata.csv \
  -c /STORAGE/csbig/sc_ADers/metadata/tables/clinical_stratified_s4.csv \
  -o /STORAGE/csbig/sc_ADers/merge_integration_results_minimal_prueba_05-22 \
  -w 4
