
python ~/Single_cell_pre/Celltypi_anno.py \
  --base /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final \
  --model /STORAGE/csbig/sc_ADers/celltipy_model_PC.pkl \
  --majority_voting \
  --min_prob 0.6 \
  --out_csv /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final/unassigned_celltypist_predictions.csv
