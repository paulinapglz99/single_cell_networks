#Plots finale for mini 

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Usage:\n  Rscript qc_plots_cowplot_style.R <RUN_DIR>\n\n",
    "Example:\n  Rscript qc_plots_cowplot_style.R \"C:/Users/danae/OneDrive/Documentos/ayuda/matrices_demultiplexed_minimal-2026-01-22_17-52\"\n"
  )
}

run_dir <- normalizePath(args[1], mustWork = TRUE)

# QC parameters FIXED in the script 
min_feat  <- 700L
min_count <- 1500L
mt_max    <- 13
down_per  <- 10000L


cat("# RUN_DIR: ", run_dir, "\n")
cat("# Thresholds: min_feat=", min_feat,
    " min_count=", min_count,
    " mt_max=", mt_max,
    " down_per=", down_per, "\n")

# Inputs from QC run
f_pre_cells <- file.path(run_dir, "qc_cells_pre.csv")
f_pre_sum   <- file.path(run_dir, "qc_pre_summary.csv")
f_post_sum  <- file.path(run_dir, "qc_post_summary.csv")
f_seq       <- file.path(run_dir, "qc_summary_cells_sequential.csv")
f_dbl_sum   <- file.path(run_dir, "doublet_summary.csv")
f_dbl_cells <- file.path(run_dir, "doublet_cells.csv")

stopifnot(file.exists(f_pre_cells))

# Output dir
plots_root <- file.path(dirname(run_dir), "qc_plots_cowplot_style")
out_dir    <- file.path(plots_root, basename(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("# Saving plots to:\n  ", out_dir, "\n")

# Save helper: PNG + JPG
save_plot <- function(plot, filename_base,
                      width = 7, height = 6, dpi = 300,
                      bg = "white") {
  ggsave(
    filename = paste0(filename_base, ".png"),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg
  )
  ggsave(
    filename = paste0(filename_base, ".jpg"),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg
  )
}

# Read PRE per-cell QC table
pre <- readr::read_csv(f_pre_cells, show_col_types = FALSE)

# Normalize percent.mt column name
if ("percent.mt" %in% names(pre) && !("percent_mt" %in% names(pre))) {
  pre <- dplyr::rename(pre, percent_mt = `percent.mt`)
}

# Ensure expected columns
need <- c("sample_id","cell_barcode","nFeature_RNA","nCount_RNA","percent_mt")
miss <- setdiff(need, names(pre))
if (length(miss)) stop("Missing in qc_cells_pre.csv: ", paste(miss, collapse=", "))

# Compute keep/remove flags if QC didn't write them
if (!("keep_cf" %in% names(pre))) {
  pre <- pre %>% mutate(keep_cf = if_else(nFeature_RNA >= min_feat & nCount_RNA >= min_count, "keep", "remove"))
}
if (!("keep_mt" %in% names(pre))) {
  pre <- pre %>% mutate(keep_mt = if_else(percent_mt <= mt_max, "keep", "remove"))
}
if (!("keep_all" %in% names(pre))) {
  pre <- pre %>% mutate(keep_all = if_else(keep_cf == "keep" & keep_mt == "keep", "keep", "remove"))
}


pre <- pre %>%
  mutate(
    keep_cf  = factor(keep_cf,  levels = c("keep","remove")),
    keep_mt  = factor(keep_mt,  levels = c("keep","remove")),
    keep_all = factor(keep_all, levels = c("keep","remove"))
  )

# Downsample per sample for plotting speed
downsample_per_sample <- function(df, by = "sample_id", cap = 10000L) {
  cap <- as.integer(cap)
  if (is.na(cap) || cap <= 0L) return(df)
  parts <- split(df, df[[by]], drop = TRUE)
  sampled <- lapply(parts, function(dd) {
    m <- min(cap, nrow(dd))
    if (m >= nrow(dd)) return(dd)
    dd[sample.int(nrow(dd), size = m), , drop = FALSE]
  })
  dplyr::bind_rows(sampled)
}
set.seed(123)
pre_ds <- downsample_per_sample(pre, by="sample_id", cap=down_per)

# Cowplot theme
base_theme <- theme_cowplot() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# global colors
qc_cols <- c(keep = "cornflowerblue", remove = "brown2")
scale_qc_keep <- function(name = "") scale_color_manual(values = qc_cols, drop = FALSE, name = name)

# colors for doublets
dbl_cols <- c(singlet = "cornflowerblue", doublet = "brown2")

# Histogram: percent_mt (log10 x-axis + ticks ~10)
p_hist_mt <- ggplot(pre, aes(x = percent_mt + 1)) +
  geom_histogram(binwidth = 0.25, fill = "cornflowerblue", colour = "black") +
  geom_vline(xintercept = mt_max + 1, linetype = 2, color = "red", linewidth = 0.8) +
  scale_x_log10(
    breaks = c(1, 11, 21, 31, 41, 51, 61),
    labels = c("0", "10", "20", "30", "40", "50", "60")
  ) +
  ggtitle("Distribution of mitochondrial expression (%)") +
  labs(x = "percent_mt (log10 axis; +1 shift)", y = "cells") +
  base_theme

save_plot(p_hist_mt, file.path(out_dir, "01_hist_percent_mt_log10_ticks10"), width=7, height=5, dpi=300)

# Histogram: nFeature_RNA + 
p_hist_nf <- ggplot(pre, aes(nFeature_RNA)) +
  geom_histogram(bins = 80, fill="cornflowerblue", colour="black") +
  geom_vline(xintercept = min_feat, linetype=2, color="red", linewidth=0.8) +
  ggtitle("Distribution of genes (nFeature_RNA)") +
  labs(x = "nFeature_RNA", y = "cells") +
  base_theme

save_plot(p_hist_nf, file.path(out_dir, "02_hist_nFeature_blue_redthreshold"), width=7, height=5, dpi=300)

# Histogram: nCount_RNA (log10) 
p_hist_nc <- ggplot(pre, aes(nCount_RNA)) +
  geom_histogram(bins = 80, fill="cornflowerblue", colour="black") +
  geom_vline(xintercept = min_count, linetype=2, color="red", linewidth=0.8) +
  scale_x_log10() +
  ggtitle("Distribution of total counts (nCount_RNA, log10)") +
  labs(x = "nCount_RNA (log10)", y = "cells") +
  base_theme

save_plot(p_hist_nc, file.path(out_dir, "03_hist_nCount_log10_blue_redthreshold"), width=7, height=5, dpi=300)

# FeatureScatter-like: nCount vs nFeature colored by percent_mt
p_mito_grad <- pre_ds %>%
  arrange(percent_mt) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour = percent_mt)) +
  geom_point(size=0.7, alpha=0.8) +
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Mitochondrial expression level") +
  labs(x = "nCount_RNA", y = "nFeature_RNA", colour = "percent_mt") +
  base_theme

save_plot(p_mito_grad, file.path(out_dir, "04_scatter_mito_gradient"), width=6.8, height=6.2, dpi=300)

# FeatureScatter-like: nCount vs nFeature colored by percent_mt (gradient, log10)
p_mito_grad_log <- pre_ds %>%
  arrange(percent_mt) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour = percent_mt)) +
  geom_point(size=0.7, alpha=0.8) +
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle("Mitochondrial expression level (log10)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)", colour = "percent_mt") +
  base_theme

save_plot(p_mito_grad_log, file.path(out_dir, "05_scatter_mito_gradient_log10"), width=6.8, height=6.2, dpi=300)

# FeatureScatter-like: nCount vs percent_mt colored by keep_mt
p_nc_mt <- ggplot(pre_ds, aes(x = nCount_RNA, y = percent_mt, color = keep_mt)) +
  geom_point(size=0.7, alpha=0.7) +
  geom_hline(yintercept = mt_max, linetype=2, color="black", linewidth=0.7) +
  scale_x_log10() +
  scale_qc_keep(name = "") +
  ggtitle("nCount_RNA vs percent_mt (log10 x) — keep/remove by MT threshold") +
  labs(x = "nCount_RNA (log10)", y = "percent_mt") +
  base_theme

save_plot(p_nc_mt, file.path(out_dir, "06_scatter_nCount_vs_percentmt_keepMT_blackthreshold"), width=7.2, height=4.8, dpi=300)

# FeatureScatter-like: nCount vs nFeature colored by keep_cf 
p_sc_keep_cf <- ggplot(pre_ds, aes(x = nCount_RNA, y = nFeature_RNA, color = keep_cf)) +
  geom_point(size=0.7, alpha=0.6) +
  geom_vline(xintercept = min_count, linetype=2, color="black", linewidth=0.6) +
  geom_hline(yintercept = min_feat,  linetype=2, color="black", linewidth=0.6) +
  scale_x_log10() +
  scale_y_log10() +
  scale_qc_keep(name = "") +
  ggtitle("Keep vs remove (counts + features thresholds)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme

save_plot(p_sc_keep_cf, file.path(out_dir, "07_scatter_keep_remove_CF_blackthresholds"), width=6.8, height=6.2, dpi=300)

# FeatureScatter-like: nCount vs nFeature colored by keep_mt
p_sc_keep_mt <- ggplot(pre_ds, aes(x = nCount_RNA, y = nFeature_RNA, color = keep_mt)) +
  geom_point(size=0.7, alpha=0.6) +
  scale_x_log10() +
  scale_y_log10() +
  scale_qc_keep(name = "") +
  ggtitle("Keep vs remove (mitochondrial threshold)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme

save_plot(p_sc_keep_mt, file.path(out_dir, "08_scatter_keep_remove_MT_colors_fixed"), width=6.8, height=6.2, dpi=300)

# FeatureScatter-like: nCount vs nFeature colored by keep_all
p_sc_keep_all <- ggplot(pre_ds, aes(x = nCount_RNA, y = nFeature_RNA, color = keep_all)) +
  geom_point(size=0.7, alpha=0.6) +
  geom_vline(xintercept = min_count, linetype=2, color="black", linewidth=0.6) +
  geom_hline(yintercept = min_feat,  linetype=2, color="black", linewidth=0.6) +
  scale_x_log10() +
  scale_y_log10() +
  scale_qc_keep(name = "") +
  ggtitle("Keep vs remove (CF + MT combined)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme

save_plot(p_sc_keep_all, file.path(out_dir, "09_scatter_keep_remove_ALL"), width=6.8, height=6.2, dpi=300)

# Violin: percent_mt by sample 
p_vln_mt <- ggplot(pre_ds, aes(x = sample_id, y = percent_mt)) +
  geom_violin(fill="grey85", colour="black", scale="width") +
  geom_boxplot(width=0.12, outlier.size=0.2, alpha=0.8) +
  geom_hline(yintercept = mt_max, linetype=2, color="red", linewidth=0.7) +
  coord_flip() +
  ggtitle("percent_mt by sample (PRE) + threshold") +
  labs(x = "sample", y = "percent_mt") +
  base_theme

save_plot(p_vln_mt, file.path(out_dir, "10_violin_percentmt_threshold_red"), width=10, height=12, dpi=300)

# Violin: nFeature_RNA by sample 
p_vln_nf <- ggplot(pre_ds, aes(x = sample_id, y = nFeature_RNA)) +
  geom_violin(fill="grey85", colour="black", scale="width") +
  geom_boxplot(width=0.12, outlier.size=0.2, alpha=0.8) +
  geom_hline(yintercept = min_feat, linetype=2, color="red", linewidth=0.7) +
  coord_flip() +
  ggtitle("nFeature_RNA by sample (PRE) + threshold") +
  labs(x = "sample", y = "nFeature_RNA") +
  base_theme

save_plot(p_vln_nf, file.path(out_dir, "11_violin_nFeature_threshold_red"), width=10, height=12, dpi=300)

# Violin: nCount_RNA by sample 
p_vln_nc <- ggplot(pre_ds, aes(x = sample_id, y = nCount_RNA)) +
  geom_violin(fill="grey85", colour="black", scale="width") +
  geom_boxplot(width=0.12, outlier.size=0.2, alpha=0.8) +
  geom_hline(yintercept = min_count, linetype=2, color="red", linewidth=0.7) +
  coord_flip() +
  ggtitle("nCount_RNA by sample (PRE) + threshold") +
  labs(x = "sample", y = "nCount_RNA") +
  base_theme

save_plot(p_vln_nc, file.path(out_dir, "12_violin_nCount_threshold_red"), width=10, height=12, dpi=300)

# Doublets: fraction per sample barplot 
if (file.exists(f_dbl_sum)) {
  dbl <- read_csv(f_dbl_sum, show_col_types = FALSE)
  if (all(c("sample","frac_doublet") %in% names(dbl))) {
    dbl <- dbl %>% mutate(frac_doublet = as.numeric(frac_doublet)) %>% arrange(desc(frac_doublet))
    p_dbl_bar <- ggplot(dbl, aes(x = reorder(sample, frac_doublet), y = frac_doublet)) +
      geom_col(fill="cornflowerblue", colour="black") +
      coord_flip() +
      scale_y_continuous(labels=percent_format(accuracy=0.1)) +
      ggtitle("Doublet fraction per sample") +
      labs(x = "sample", y = "doublet fraction") +
      base_theme
    save_plot(p_dbl_bar, file.path(out_dir, "13_doublet_fraction_per_sample"), width=8, height=10, dpi=300)
  }
}

# FeatureScatter-like: nCount vs nFeature colored by scDblFinder.class 
if (file.exists(f_dbl_cells)) {
  dbl_cells <- read_csv(f_dbl_cells, show_col_types = FALSE)
  
  if ("barcode" %in% names(dbl_cells) && !("cell_barcode" %in% names(dbl_cells))) {
    dbl_cells <- rename(dbl_cells, cell_barcode = barcode)
  }
  
  dbl_plot <- pre %>%
    select(sample_id, cell_barcode, nCount_RNA, nFeature_RNA) %>%
    left_join(
      dbl_cells %>% select(sample_id, cell_barcode, scDblFinder.class),
      by = c("sample_id","cell_barcode")
    ) %>%
    filter(!is.na(scDblFinder.class)) %>%
    mutate(
      scDblFinder.class = factor(scDblFinder.class, levels = c("singlet","doublet")),
      keep_dbl = factor(if_else(scDblFinder.class == "doublet", "remove", "keep"),
                        levels = c("keep","remove"))
    )
  
  set.seed(123)
  dbl_plot_ds <- downsample_per_sample(dbl_plot, by="sample_id", cap=down_per)
  
  #  keep/remove (doublets)
  p_sc_doublets_keepremove <- ggplot(dbl_plot_ds,
                                     aes(x = nCount_RNA, y = nFeature_RNA, color = keep_dbl)) +
    geom_point(size=0.7, alpha=0.6) +
    geom_vline(xintercept = min_count, linetype=2, color="black", linewidth=0.6) +
    geom_hline(yintercept = min_feat,  linetype=2, color="black", linewidth=0.6) +
    scale_x_log10() +
    scale_y_log10() +
    scale_qc_keep(name = "") +
    ggtitle("Doublets QC — keep vs remove (scDblFinder)") +
    labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
    base_theme
  
  save_plot(p_sc_doublets_keepremove, file.path(out_dir, "14_scatter_doublets_keep_remove_colors_fixed"), width=6.8, height=6.2, dpi=300)
  
  # singlet vs doublet 
  p_sc_doublets_class <- ggplot(dbl_plot_ds,
                                aes(x = nCount_RNA, y = nFeature_RNA, color = scDblFinder.class)) +
    geom_point(size=0.7, alpha=0.6) +
    geom_vline(xintercept = min_count, linetype=2, color="black", linewidth=0.6) +
    geom_hline(yintercept = min_feat,  linetype=2, color="black", linewidth=0.6) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = dbl_cols, drop = FALSE, name = "") +
    ggtitle("scDblFinder classification — singlets vs doublets") +
    labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
    base_theme
  
  save_plot(p_sc_doublets_class, file.path(out_dir, "15_scatter_doublets_singlet_vs_doublet_colors_fixed"), width=6.8, height=6.2, dpi=300)
}

# PRE vs POST: merge summaries 
if (file.exists(f_pre_sum) && file.exists(f_post_sum)) {
  preS  <- read_csv(f_pre_sum, show_col_types = FALSE)
  postS <- read_csv(f_post_sum, show_col_types = FALSE)
  
  pp <- full_join(preS, postS, by="sample", suffix=c("_pre","_post"))
  write_csv(pp, file.path(out_dir, "qc_pre_post_merged.csv"))
  
  if (all(c("cells_pre","cells_post") %in% names(pp))) {
    p_cells <- pp %>%
      pivot_longer(c(cells_pre, cells_post), names_to="stage", values_to="cells") %>%
      ggplot(aes(x=reorder(sample, cells), y=cells, fill=stage)) +
      geom_col(position="dodge", colour="black") +
      coord_flip() +
      ggtitle("Cells per sample (PRE vs POST)") +
      labs(x="sample", y="cells") +
      base_theme
    save_plot(p_cells, file.path(out_dir, "16_cells_pre_vs_post"), width=9, height=12, dpi=300)
  }
}

# 2 plots per page
add_if_exists <- function(lst, obj_name) {
  if (exists(obj_name, inherits = TRUE)) {
    lst[[length(lst) + 1]] <- get(obj_name, inherits = TRUE)
  }
  lst
}

plots_to_pdf <- list()

plots_to_pdf <- add_if_exists(plots_to_pdf, "p_hist_mt")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_hist_nf")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_hist_nc")

plots_to_pdf <- add_if_exists(plots_to_pdf, "p_mito_grad")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_mito_grad_log")

plots_to_pdf <- add_if_exists(plots_to_pdf, "p_nc_mt")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_sc_keep_cf")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_sc_keep_mt")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_sc_keep_all")

plots_to_pdf <- add_if_exists(plots_to_pdf, "p_vln_mt")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_vln_nf")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_vln_nc")

plots_to_pdf <- add_if_exists(plots_to_pdf, "p_sc_doublets_keepremove")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_sc_doublets_class")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_dbl_bar")
plots_to_pdf <- add_if_exists(plots_to_pdf, "p_cells")

if (length(plots_to_pdf) == 0) {
  stop("No plots found to write to PDF. (plots_to_pdf is empty)")
}

pdf_file <- file.path(out_dir, "QC_all_plots_2perpage.pdf")

pdf(pdf_file, width = 8.5, height = 11)  
for (i in seq(1, length(plots_to_pdf), by = 2)) {
  p1 <- plots_to_pdf[[i]]
  p2 <- if (i + 1 <= length(plots_to_pdf)) plots_to_pdf[[i + 1]] else NULL
  
  if (!is.null(p2)) {
    page <- cowplot::plot_grid(p1, p2, nrow = 2, ncol = 1, align = "v")
  } else {
    page <- cowplot::plot_grid(p1, nrow = 1, ncol = 1)
  }
  print(page)
}
dev.off()

cat("# Saved PDF: ", pdf_file, "\n")
cat("# Total plots written: ", length(plots_to_pdf), "\n")
cat("# Done. Plots in:\n  ", out_dir, "\n")
