<!-- ===================== HEADER (50/50) ===================== -->
<table>
  <tr>
    <td width="55%" valign="top">

<h1>Single cell networks repository</h1>

<p>
This repository implements a reproducible pipeline for single-cell RNA sequencing (scRNA-seq) analysis, following current best practices for
single-cell preprocessing as described in Luecken &amp; Theis (2019).
</p>

<p><b>The primary objectives of this project are:</b></p>

<ul>
  <li>
    To perform standardized preprocessing of single-cell data, including metadata harmonization, demultiplexing, quality control, normalization, integration, and cell annotation.
  </li>
  <li>
    To construct gene co-expression networks from processed single-cell data for downstream biological analysis.
  </li>
</ul>

<p>
<b>Language:</b> R<br/>
<b>Core library:</b> Seurat<br/>
<b>Additional packages include:</b> Harmony, scDblFinder, tidyverse, SingleCellExperiment, etc
</p>

    </td>
    <td width="45%" align="center" valign="top">

<img width="40%" alt="Copy of worflow_for_git drawio (2)"
     src="https://github.com/user-attachments/assets/4bc97aa7-22fe-488c-8b0b-4145c2861fa6" />

    </td>
  </tr>
</table>

<hr/>

<!-- ===================== DATA ACQUISITION ===================== -->
<h2>Data Acquisition</h2>

<p>All data were obtained from Synapse with approved access.</p>

<p>
<b>Study:</b> ROSMAP – DLPFC Experiment 2<br/>
<b>Two types of data were available on Synapse:</b>
</p>

<ul>
  <li>Raw FASTQ files</li>
  <li>Processed count matrices - <code>syn3157322</code></li>
</ul>

<h3>Metadata Files</h3>

<p>The following metadata files were downloaded from Synapse: <code>syn3157322</code></p>

<ul>
  <li>ROSMAP_clinical.csv</li>
  <li>ROSMAP_biospecimen_metadata.csv</li>
  <li>ROSMAP_assay_scrnaSeq_metadata.csv</li>
  <li>ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv</li>
</ul>

<hr/>

<!-- ===================== DEMULTIPLEXING ===================== -->
<h2>Demultiplexing</h2>

<p><b>Script :</b> <code>1.demultiplex_matrices.R</code></p>

<p>
<b>Background: why demultiplexing is required in ROSMAP (DLPFC Experiment 2) :</b>
In this dataset, most count matrices correspond to pooled libraries, meaning that a single matrix contains cells from multiple donors.
</p>

<p>Libraries were generated from pooled samples as follows:</p>

<ul>
  <li>222 libraries include 8 donors</li>
  <li>4 libraries include 7 donors</li>
  <li>8 libraries include a single donor</li>
</ul>

<p><b>Additionally:</b></p>
<ul>
  <li>Most pooled libraries were prepared as two replicate library preparations (e.g., B10-A and B10-B).</li>
  <li>Each replicate library was sequenced at two different sequencing centers (Broad and NYGC).</li>
</ul>

<p>Therefore, a library batch (e.g., B10) can produce four sequencing datasets:</p>

<ul>
  <li>B10-A-Broad, B10-A-NYGC, B10-B-Broad, and B10-B-NYGC.</li>
</ul>

<p>Because each pooled matrix contains multiple donors, we must assign each cell barcode to its donor before any other process.</p>

<table>
  <tr>
    <td width="50%" valign="top">
      <h3>Input</h3>
      <ul>
        <li>Input : count matrices (<code>syn3157322</code>)</li>
      </ul>

      <p>To do this, we use the demultiplexing mapping file:</p>
      <p><code>ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv (syn34572333)</code></p>

      <p>This file links:</p>
      <p><code>cellBarcode + libraryBatch → individualID</code></p>

      <ul>
        <li>Input : count matrices (<code>syn3157322</code>)</li>
      </ul>
    </td>

    <td width="50%" valign="top">
      <h3>Output</h3>
      <p>
        Output file : <code>matrices_demultiplexed_final.rds</code>
        → A flattened list of Seurat objects, where each element corresponds to an individual donor within a given library batch (named as <code>libraryBatch_individualID</code>).
      </p>

      <h3>How to run</h3>
      <p><b>Command :</b> <code>bash run_1.demultiplex.sh</code></p>

      <p><code>Rscript .../CopyOf1.demultiplex_matrices.R</code> -> Runs the demultiplexing R script</p>

      <ul>
        <li><code>-directory /datos/rosmap/single_cell/matrix_exp_2/</code> -> Path to the folder containing the input count matrices downloaded from Synapse (Experiment 2).</li>
        <li><code>--metadata .../ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv</code> -> Path to the demultiplexing mapping file. This file provides the key mapping: cellBarcode + libraryBatch → individualID.</li>
        <li><code>--output .../matrices_demultiplexed_proof.rds</code> -> Output path for the generated .rds object.</li>
        <li><code>--test</code> -> Runs the script in test mode (typically used to run a smaller subset / quick validation).</li>
      </ul>
    </td>
  </tr>
</table>

<hr/>

<!-- ===================== QUALITY CONTROL ===================== -->
<h2>Quality control</h2>

<p><b>Script :</b> <code>2.quality_control.R</code></p>

<p>
Make structural integrity validation of the input Seurat objects (Seurat v5 layers-aware), ensuring the data are readable, consistent, and properly annotated before filtering.
</p>

<p>
Perform standardized, multi-metric quality control on demultiplexed Seurat objects to remove low-quality cells, lowly detected genes, and predicted doublets before normalization/integration.
</p>

<table>
  <tr>
    <td width="50%" valign="top">
      <h3>Input</h3>
      <ul>
        <li>Input - A list of Seurat objects -> <code>matrices_demultiplexed_final.rds</code></li>
      </ul>
    </td>
    <td width="50%" valign="top">
      <h3>Output</h3>
      <p>A timestamped output folder created in the current working directory:</p>
      <pre><code>&lt;INPUT_BASENAME&gt;-YYYY-MM-DD_HH-MM/</code></pre>

      <p><b>Main outputs include:</b></p>
      <ul>
        <li><code>seurat_list_filtered.rds</code> — filtered Seurat object list (ready for normalization/integration)</li>
        <li><code>damaged_samples.txt</code> — samples dropped due to structural issues (invalid assay layers or missing QC metrics)</li>
        <li><code>qc_cells_pre.csv</code> — per-cell QC metrics before filtering</li>
        <li><code>qc_pre_summary.csv</code> — per-sample summary statistics before filtering</li>
        <li><code>doublet_cells.csv</code> — per-cell scDblFinder calls</li>
        <li><code>doublet_summary.csv</code> — per-sample doublet rate summary</li>
        <li><code>qc_post_summary.csv</code> — per-sample summary after QC filtering</li>
        <li><code>qc_summary_cells_sequential.csv</code> — total cell counts retained after each QC step</li>
      </ul>
    </td>
  </tr>
</table>

<p>
<b>QC thresholds (current defaults):</b>
Cells are retained if they have at least 700 detected genes (nFeature_RNA ≥ 700) and 1500 UMIs (nCount_RNA ≥ 1500), with a maximum mitochondrial content of 13% (percent.mt ≤ 13).
Genes are kept only if detected in at least 10 cells across the dataset.
Doublet detection is enabled by default using scDblFinder, and mitochondrial genes are identified using the human gene prefix pattern "^MT-".
</p>

<h3>How to run</h3>
<ul>
  <li>Command : <code>bash run_2.quality_control.sh</code></li>
</ul>

<p><code>Rscript ~/single_cell_networks/2.quality_control.R</code> -> Runs the QC script.</p>
<p><code>/STORAGE/csbig/sc_ADers/matrices_demultiplexed_final.rds</code> -> (input): path to the demultiplexed Seurat list produced in demultiplexing step.</p>

<hr/>

<!-- ===================== QC PLOTS ===================== -->
<h2>QC Plots</h2>

<p><b>Script :</b> <code>2.1.qc_plots.R</code></p>

<p>
This script takes as input the QC run output directory (the folder created by 2.quality_control.R) and reads the summary tables generated during QC .
Using these files, it produces a complete set of QC plots to evaluate filtering thresholds and sample-level quality.
</p>

<table>
  <tr>
    <td width="50%" valign="top">
      <h3>Input</h3>
      <p>A QC output folder created in Step 2, for example:</p>
      <pre><code>matrices_demultiplexed_final-YYYY-MM-DD_HH-MM/</code></pre>
    </td>

    <td width="50%" valign="top">
      <h3>Output</h3>
      <p>A plots directory automatically created at:</p>
      <pre><code>qc_plots_cowplot_style/&lt;QC_RUN_FOLDER_NAME&gt;/</code></pre>

      <p>This folder contains:</p>
      <ul>
        <li>Multiple QC plots saved in PNG and JPG</li>
      </ul>

      <p>A combined PDF report:</p>
      <ul>
        <li><code>QC_all_plots_2perpage.pdf</code></li>
      </ul>
    </td>
  </tr>
</table>

<h3>How to run</h3>
<ul>
  <li>Command : <code>bash run_2.1_qc_plots.sh</code></li>
</ul>

<p><code>Rscript ~/single_cell_networks/2.1.qc_plots.R</code> -> Runs the plotting script.</p>

<p>
<code>/STORAGE/.../matrices_demultiplexed_final_QC-2026-01-28_19-55</code>
-> (input): the QC output folder generated in Step 2. The script expects QC summary tables inside this folder and will create plots accordingly.
</p>

<hr/>

<!-- ===================== NORMALIZATION / MERGE / INTEGRATION ===================== -->
<h2>Normalization, Merge, and Integration</h2>

<p><b>Script :</b> <code>3.merge_integration.R</code></p>

<p>
Merge QC-filtered Seurat objects, perform standard log-normalization, and integrate the dataset while correcting batch effects using Harmony.
This step generates an integrated embedding (UMAP), clustering results, and a final merged Seurat object
</p>

<table>
  <tr>
    <td width="50%" valign="top">
      <h3>Input</h3>
      <p><b>Input :</b> <code>seurat_list_filtered.rds</code> (output from Step 2: QC-filtered Seurat object list)</p>
      <p>ROSMAP_assay_scrnaSeq_metadata.csv, clinical_stratified.csv</p>
    </td>

    <td width="50%" valign="top">
      <h3>Output</h3>
      <p>
        An output directory (set by <code>--out_dir</code>) containing:
      </p>
      <ul>
        <li><code>merged_by_individual_harmony.rds</code> : Final integrated Seurat object (Harmony-corrected)</li>
        <li><code>merged_by_individual_harmony_cell_metadata.csv</code>, Cell-level metadata from the final object</li>
        <li>UMAPs and graphics</li>
      </ul>
    </td>
  </tr>
</table>

<h3>Notes</h3>
<ul>
  <li>Normalization method: Seurat NormalizeData() (log-normalization).</li>
  <li>Batch correction: Harmony integration on libraryBatch (other variables like platformLocation can be added if needed).</li>
</ul>

<h3>How to run</h3>
<ul>
  <li>Command : <code>bash run_3.merge_integration.sh</code></li>
</ul>

<p><code>Rscript ~/single_cell_networks/3.merge_integration.R</code> -> Runs the merge + normalization + Harmony integration script.</p>

<ul>
  <li><code>-s .../seurat_list_filtered.rds</code> -> (input) : QC-filtered Seurat list produced by qc , This is the main object list that will be merged and integrated.</li>
  <li><code>-a .../ROSMAP_assay_scrnaSeq_metadata.csv</code> -> assay-level metadata. Used to attach/validate experiment-level variables (e.g., library batch, sequencing batch, center/platform).</li>
  <li><code>-c .../clinical_stratified.csv</code> -> Clinical metadata table (preprocessed/stratified). Used to append donor-level clinical covariates (e.g., diagnosis groups, demographics, etc.) to the merged object.</li>
  <li><code>-o .../merge_integration_results_feb10</code> -> Output directory where all integration results will be saved (final .rds, plots, exported metadata).</li>
  <li><code>-w 12</code> -> Number of workers/threads used for parallel steps</li>
</ul>











