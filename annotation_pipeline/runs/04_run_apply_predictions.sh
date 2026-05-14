Rscript ~/single_cell_networks/04_apply_predictions_CT.R\
   --seurat_rds /STORAGE/csbig/sc_ADers/out_annotated_final_s4/merged_harmony_integrated_annotated.rds \
   --pred_csv /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final_s4/unassigned_celltypist_predictions_s4.csv \
   --out_dir /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final_s4/annotation_with_celltypist_s4 \
   -w 12 \
   --seed 42
