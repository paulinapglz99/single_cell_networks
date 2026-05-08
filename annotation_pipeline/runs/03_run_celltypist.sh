
python ~/single_cell_networks/03_Celltypist_anno.py \
  --base /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final_s4 \
  --model /STORAGE/csbig/sc_ADers/celltipy_model_PC.pkl \
  --majority_voting \
  --min_prob 0.6 \
  --out_csv /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_final_s4/unassigned_celltypist_predictions_s4.csv
