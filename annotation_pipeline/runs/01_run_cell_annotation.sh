Rscript ~/Single_cell_pre/cell_annotation.R \
  --seurat_rds /STORAGE/csbig/sc_ADers/merge_integration_results_feb10/merged_by_individual_harmony.rds \
  --anno_csv /STORAGE/csbig/sc_ADers/metadata/Experiment2/cell-annotation.full-atlas.csv \
  --out_dir /STORAGE/csbig/sc_ADers/out_annotated_final \
  -w 12 \
  --seed 42
