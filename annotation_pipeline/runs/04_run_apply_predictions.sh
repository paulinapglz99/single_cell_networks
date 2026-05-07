Rscript ~/Single_cell_pre/apply_predictions_CT.R \
   --seurat_rds /STORAGE/csbig/sc_ADers/out_annotated_final/merged_harmony_integrated_annotated.rds \
   --pred_csv /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final/unassigned_celltypist_predictions.csv \
   --out_dir /STORAGE/csbig/sc_ADers/celltypist/celltypist/unassigned_only_final/annotation_with_celltypist \
   -w 12 \
   --seed 42
