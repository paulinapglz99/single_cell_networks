#!/bin/bash
Rscript ~/single_cell_networks/6_builmatrices_runbash.R \
  --seurat_rds /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_minimal_2/annotation_with_celltypist/merged_harmony_integrated_annotated_plus_celltypist.rds \
  --sc_rds     /STORAGE/csbig/sc_ADers/supercell/supercell_minimal_26_05/SC_gamma20_2026-05-26.rds \
  --out_dir    /STORAGE/csbig/sc_ADers/supercell/supercell_minimal_26_05/matrices_aracne_try2 \
  --cell_types "Excitatory Neurons,Inhibitory Neurons,Astrocyte,Oligodendrocytes,Endothelial,OPCs,Microglia" \
  --slot       data
