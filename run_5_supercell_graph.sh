Rscript ~/single_cell_networks/5.supercell_graph.R \
  --seurat_rds /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_minimal_2/annotation_with_celltypist/merged_harmony_integrated_annotated_plus_celltypist.rds \
  --out_dir    /STORAGE/csbig/sc_ADers/supercell_test \
  --gammas     "10,20,50" \
  --cell_type_col cell_type \
  --phenotype_col is_AD \
  --reduction  harmony \
  --n_pc       30 \
  --k_knn      5 \
  --seed       42
