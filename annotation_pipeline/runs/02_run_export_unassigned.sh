Rscript ~/Single_cell_pre/select_na.R \
  --seurat_rds /STORAGE/csbig/sc_ADers/out_annotated_final/merged_harmony_integrated_annotated.rds \
  --out_dir /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final \
  -w 12 \
  --seed 42
