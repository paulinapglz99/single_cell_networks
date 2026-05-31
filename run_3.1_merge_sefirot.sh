#!/bin/bash

Rscript ~/single_cell_networks/3.merge_joinlayers \
  -s /datos/rosmap/single_cell/data_from_fenix2026-05-28/matrices_demultiplexed_final_QC-2026-01-28_19-55/seurat_list_filtered.rds \
  -a /datos/rosmap/single_cell/data_from_fenix2026-05-28/metadata/ROSMAP_assay_scrnaSeq_metadata.csv \
  -c /datos/rosmap/single_cell/data_from_fenix2026-05-28/metadata/tables/clinical_stratified_s4.csv \
  -o /datos/rosmap/single_cell/merge_integration_joinlayers_2026-05-31 \
  -w 2
