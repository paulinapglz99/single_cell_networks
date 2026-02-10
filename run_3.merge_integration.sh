Rscript ~/single_cell_networks/3.merge_integration.R \
-s /STORAGE/csbig/sc_ADers/matrices_demultiplexed_final_QC-2026-01-28_19-55/seurat_list_filtered.rds \
-a /STORAGE/csbig/sc_ADers/metadata/ROSMAP_assay_scrnaSeq_metadata.csv \
-c /STORAGE/csbig/sc_ADers/metadata/tables/clinical_stratified.csv \
-o /STORAGE/csbig/sc_ADers/merge_integration_results_feb10 \
-w 12
