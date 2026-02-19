# Single cell networks repository 

This repository implements a reproducible pipeline for single-cell RNA sequencing (scRNA-seq) analysis, following current best practices for single-cell preprocessing as described in Luecken & Theis (2019).
The primary objectives of this project are:
To perform standardized preprocessing of single-cell data, including metadata harmonization, demultiplexing, quality control, normalization, integration, and cell annotation.


To construct gene co-expression networks from processed single-cell data for downstream biological analysis.


Language: R


Core library: Seurat


Additional packages include: Harmony, scDblFinder, tidyverse, SingleCellExperiment, etc 

 ## Data Acquisition
All data were obtained from Synapse with approved access.
Study: ROSMAP – DLPFC Experiment 2
Two types of data were available on Synapse:
Raw FASTQ files


Processed count matrices - syn3157322



Metadata Files
The following metadata files were downloaded from Synapse: syn3157322

ROSMAP_clinical.csv


ROSMAP_biospecimen_metadata.csv


ROSMAP_assay_scrnaSeq_metadata.csv


ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv


## Demultiplexing 
Script : 1.demultiplex_matrices.R

Background: why demultiplexing is required in ROSMAP (DLPFC Experiment 2) : In this dataset, most count matrices correspond to pooled libraries, meaning that a single matrix contains cells from multiple donors.
 Libraries were generated from pooled samples as follows: 222 libraries include 8 donors, 4 libraries include 7 donors, and 8 libraries include a single donor.
Additionally:
Most pooled libraries were prepared as two replicate library preparations (e.g., B10-A and B10-B).


Each replicate library was sequenced at two different sequencing centers (Broad and NYGC).


Therefore, a library batch (e.g., B10) can produce four sequencing datasets:
 B10-A-Broad, B10-A-NYGC, B10-B-Broad, and B10-B-NYGC.


Because each pooled matrix contains multiple donors, we must assign each cell barcode to its donor before any other process.

Input : count matrices (syn3157322) 


To do this, we use the demultiplexing mapping file:
ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv (syn34572333)


This file links:
cellBarcode + libraryBatch → individualID

- Input : count matrices (syn3157322) 

- Output file  : matrices_demultiplexed_final.rds  -> A flattened list of Seurat objects, where each element corresponds to an individual donor within a given library batch (named as libraryBatch_individualID).

### How to run 
- Command : bash run_1.demultiplex.sh
  
Rscript .../CopyOf1.demultiplex_matrices.R -> Runs the demultiplexing R script

-directory /datos/rosmap/single_cell/matrix_exp_2/ -> Path to the folder containing the input count matrices downloaded from Synapse (Experiment 2).

--metadata .../ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv -> Path to the demultiplexing mapping file. This file provides the key mapping: cellBarcode + libraryBatch → individualID.

--output .../matrices_demultiplexed_proof.rds -> Output path for the generated .rds object.

--test -> Runs the script in test mode (typically used to run a smaller subset / quick validation).


## Quality control 

Script : 2.quality_control.R 

Make  structural integrity validation of the input Seurat objects (Seurat v5 layers-aware), ensuring the data are readable, consistent, and properly annotated before filtering.


Perform standardized, multi-metric quality control on demultiplexed Seurat objects to remove low-quality cells, lowly detected genes, and predicted doublets before normalization/integration.

- Input - A list of Seurat objects -> matrices_demultiplexed_final.rds

- Output
A timestamped output folder created in the current working directory:

<INPUT_BASENAME>-YYYY-MM-DD_HH-MM/

Main outputs include:
seurat_list_filtered.rds — filtered Seurat object list (ready for normalization/integration)
damaged_samples.txt — samples dropped due to structural issues (invalid assay layers or missing QC metrics)
qc_cells_pre.csv — per-cell QC metrics before filtering
qc_pre_summary.csv — per-sample summary statistics before filtering
doublet_cells.csv — per-cell scDblFinder calls 
doublet_summary.csv — per-sample doublet rate summary 
qc_post_summary.csv — per-sample summary after QC filtering
qc_summary_cells_sequential.csv — total cell counts retained after each QC step


QC thresholds (current defaults): Cells are retained if they have at least 700 detected genes (nFeature_RNA ≥ 700) and 1500 UMIs (nCount_RNA ≥ 1500), with a maximum mitochondrial content of 13% (percent.mt ≤ 13). Genes are kept only if detected in at least 10 cells across the dataset. Doublet detection is enabled by default using scDblFinder, and mitochondrial genes are identified using the human gene prefix pattern "^MT-".

### How to run 
Command : bash run_2.quality_control.sh

Rscript ~/single_cell_networks/2.quality_control.R -> Runs the QC script.

/STORAGE/csbig/sc_ADers/matrices_demultiplexed_final.rds -> (input): path to the demultiplexed Seurat list produced in demultiplexing step. 



## QC Plots
Script : 2.1.qc_plots.R

This script takes as input the QC run output directory (the folder created by 2.quality_control.R) and reads the summary tables generated during QC . Using these files, it produces a complete set of QC plots to evaluate filtering thresholds and sample-level quality.

- Input
A QC output folder created in Step 2, for example:
matrices_demultiplexed_final-YYYY-MM-DD_HH-MM/
- Output
A plots directory automatically created at:
qc_plots_cowplot_style/<QC_RUN_FOLDER_NAME>/

This folder contains:
Multiple QC plots saved in PNG and JPG

A combined PDF report:
 QC_all_plots_2perpage.pdf

### How to run 
Command : bash run_2.1_qc_plots.sh

Rscript ~/single_cell_networks/2.1.qc_plots.R -> Runs the plotting script.

/STORAGE/.../matrices_demultiplexed_final_QC-2026-01-28_19-55  -> (input): the QC output folder generated in Step 2. The script expects QC summary tables inside this folder and will create plots accordingly.

 

## Normalization, Merge, and Integration

Script : 3.merge_integration.R

Merge QC-filtered Seurat objects, perform standard log-normalization, and integrate the dataset while correcting batch effects using Harmony. This step generates an integrated embedding (UMAP), clustering results, and a final merged Seurat object 

Input :seurat_list_filtered.rds (output from Step 2: QC-filtered Seurat object list)

ROSMAP_assay_scrnaSeq_metadata.csv, clinical_stratified.csv 
- Output -> An output directory (set by --out_dir) containing:,  merged_by_individual_harmony.rds :  Final integrated Seurat object (Harmony-corrected) , merged_by_individual_harmony_cell_metadata.csv, Cell-level metadata from the final object, UMAPs and graphics

Notes
Normalization method: Seurat NormalizeData() (log-normalization).


Batch correction: Harmony integration on libraryBatch (other variables like platformLocation can be added if needed).

### How to run 
Command : bash run_3.merge_integration.sh

Rscript ~/single_cell_networks/3.merge_integration.R -> Runs the merge + normalization + Harmony integration script.

-s .../seurat_list_filtered.rds -> (input) : QC-filtered Seurat list produced by qc , This is the main object list that will be merged and integrated. 

-a .../ROSMAP_assay_scrnaSeq_metadata.csv ->  assay-level metadata. Used to attach/validate experiment-level variables (e.g., library batch, sequencing batch, center/platform).

-c .../clinical_stratified.csv ->  Clinical metadata table (preprocessed/stratified). Used to append donor-level clinical covariates (e.g., diagnosis groups, demographics, etc.) to the merged object.








