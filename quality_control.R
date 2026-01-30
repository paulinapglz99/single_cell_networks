#QC with all metrics 

#Args

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage:\n  Rscript qc_run.R <INPUT_RDS>\n\nExample:\n  Rscript qc_run.R /path/to/matrices_demultiplexed_minimal.rds\n")
}

input_path <- args[1]
stopifnot(file.exists(input_path))

# Libraries 
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(jsonlite)
  library(readr)
})

msg <- function(...) {
  cat(paste0(...), "\n")
}


msg("# Packages loaded")

#Parameters 
min_feat       <- 700L     # min nFeature_RNA per cell 
min_count      <- 1500L    # min nCount_RNA per cell  
mt_max         <- 13       # max percent.mt per cell  
min_gene_cells <- 10L      # min cells across all samples to keep a gene
run_doublets   <- TRUE     # run scDblFinder
mt_pattern     <- "^MT-"   # human: "^MT-" ; mouse: "^mt-"

# Output folder
stamp  <- format(Sys.time(), "%Y-%m-%d_%H-%M")
base   <- tools::file_path_sans_ext(basename(input_path))

base_outdir <- getwd()
outdir <- file.path(base_outdir, paste0(base, "-", stamp))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
msg("# Outdir (PWD): ", normalizePath(outdir, winslash = "/", mustWork = FALSE))

msg("# QC run started")
msg("# Input: ", input_path)
msg("# Outdir: ", outdir)
msg("# Thresholds: min_feat=", min_feat, " min_count=", min_count,
    " mt_max=", mt_max, " min_gene_cells=", min_gene_cells)
msg("# Doublets: ", run_doublets)
msg("# MT pattern: ", mt_pattern)

# Helpers
get_counts <- function(obj, assay = "RNA") {
  # Robustly get counts (Seurat v5 layers-friendly)
  if (assay %in% names(obj@assays)) {
    a <- obj[[assay]]
    if (!is.null(a@layers) && "counts" %in% names(a@layers)) {
      return(GetAssayData(obj, assay = assay, layer = "counts"))
    } else {
      return(GetAssayData(obj, assay = assay, slot = "counts"))
    }
  } else stop("Assay not found: ", assay)
}

check_assay_layers <- function(obj, assay = "RNA") {
  # Ensure assay layers are proper 2D matrices
  if (!assay %in% names(obj@assays)) return(FALSE)
  a <- obj[[assay]]
  ln <- names(a@layers)
  if (length(ln) == 0) return(FALSE)
  for (ly in ln) {
    L <- a@layers[[ly]]
    if (!(is.matrix(L) || inherits(L, "Matrix"))) return(FALSE)
    if (length(dim(L)) != 2) return(FALSE)
  }
  TRUE
}

sum_presence_across <- function(seurat_list) {
  # Sum gene presence (>0 counts) across samples without huge cbind
  total <- integer(0); names(total) <- character(0)
  for (nm in names(seurat_list)) {
    m <- get_counts(seurat_list[[nm]])
    v <- Matrix::rowSums(m > 0)
    v <- as.integer(v); names(v) <- rownames(m)
    if (length(total) == 0) {
      total <- v
    } else {
      common <- intersect(names(total), names(v))
      if (length(common) > 0) total[common] <- total[common] + v[common]
      newg <- setdiff(names(v), names(total))
      if (length(newg) > 0) total[newg] <- v[newg]
    }
  }
  total
}

msg("# Helpers ready")

# Load input and normalize to list 
x <- readRDS(input_path)
msg("# x class: ", paste(class(x), collapse = ", "))

if (inherits(x, "Seurat")) {
  seurat_list <- list(x)
  nm <- if ("sample_id" %in% colnames(x@meta.data)) unique(x$sample_id)[1] else "sample1"
  names(seurat_list) <- nm
} else if (is.list(x) && all(vapply(x, function(z) inherits(z, "Seurat"), logical(1)))) {
  seurat_list <- x
  if (is.null(names(seurat_list))) names(seurat_list) <- paste0("sample", seq_along(seurat_list))
} else {
  stop("Input .rds must be a Seurat object or a list of Seurat objects.")
}

msg("# Samples loaded: ", length(seurat_list))
msg("# Example samples: ", paste(head(names(seurat_list), 5), collapse = ", "))
msg("# Cells in first sample (raw): ", ncol(seurat_list[[1]]))

#Validate layers + compute percent.mt + ensure sample_id
bad <- character(0); i <- 0
for (nm in names(seurat_list)) {
  i <- i + 1
  if (i %% 25 == 0) msg("... validated samples: ", i)
  obj <- seurat_list[[nm]]
  
  # layers check
  if (!check_assay_layers(obj, "RNA")) { bad <- c(bad, nm); next }
  
  # compute percent.mt (we WILL filter by mt at the END)
  DefaultAssay(obj) <- "RNA"
  ok <- TRUE
  obj <- tryCatch({ obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern); obj },
                  error = function(e) { ok <<- FALSE; obj })
  if (!ok) { bad <- c(bad, nm); next }
  
  # ensure sample_id
  md <- obj@meta.data
  if (!"sample_id" %in% colnames(md)) {
    if ("sample" %in% colnames(md))          obj$sample_id <- obj$sample
    else if ("orig.ident" %in% colnames(md)) obj$sample_id <- obj$orig.ident
    else                                     obj$sample_id <- rep(nm, ncol(obj))
  }
  seurat_list[[nm]] <- obj
}

# drop damaged
if (length(bad) > 0) {
  msg("# Dropping damaged samples: ", length(bad))
  writeLines(bad, file.path(outdir, "damaged_samples.txt"))
  seurat_list <- seurat_list[setdiff(names(seurat_list), bad)]
} else msg("# No damaged samples detected")
stopifnot(length(seurat_list) > 0)

msg("# Remaining samples: ", length(seurat_list))
msg("# Cells in first sample (post-validate): ", ncol(seurat_list[[1]]))
msg("# Mean percent.mt (first sample): ", round(mean(seurat_list[[1]]$percent.mt, na.rm = TRUE), 2))
msg("# RNA layers in first sample: ", paste(names(seurat_list[[1]][["RNA"]]@layers), collapse = ", "))
msg("# Head meta.data of first sample:")
print(head(seurat_list[[1]]@meta.data[, c("sample_id","nFeature_RNA","nCount_RNA","percent.mt"), drop=FALSE]))

# Optional: drop NA in percent.mt (just in case) 
seurat_list <- lapply(seurat_list, function(obj) {
  md <- obj@meta.data
  keep <- !is.na(md$percent.mt)
  if (any(!keep)) obj <- subset(obj, cells = rownames(md)[keep])
  obj
})
msg("# NA percent.mt removed (if any)")

#  Gene filtering by global presence 
msg("# Gene filtering: computing global presence (>= min_gene_cells)")
gene_presence <- sum_presence_across(seurat_list)
msg("# Genes total in presence vector: ", length(gene_presence))
msg("# Genes with presence >= ", min_gene_cells, ": ",
    sum(gene_presence >= min_gene_cells, na.rm = TRUE))
genes_keep <- names(gene_presence[gene_presence >= min_gene_cells])
stopifnot(length(genes_keep) > 0)
seurat_list <- lapply(seurat_list, function(obj) subset(obj, features = genes_keep))
msg("# Genes kept: ", length(genes_keep))

#  Export per-cell QC table (PRE: before any cell filtering)
msg("# Exporting per-cell QC table (pre-filter)")

qc_cells_pre <- purrr::map2_dfr(
  seurat_list, names(seurat_list),
  function(obj, nm) {
    tibble::tibble(
      sample_id    = nm,
      cell_barcode = colnames(obj),
      nFeature_RNA = obj$nFeature_RNA,
      nCount_RNA   = obj$nCount_RNA,
      percent_mt   = obj$percent.mt
    )
  }
)

readr::write_csv(qc_cells_pre, file.path(outdir, "qc_cells_pre.csv"))

qc_pre_summary <- qc_cells_pre %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    cells   = dplyr::n(),
    mean_mt = mean(percent_mt,   na.rm = TRUE),
    mean_nF = mean(nFeature_RNA, na.rm = TRUE),
    mean_nC = mean(nCount_RNA,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(sample = sample_id) %>%
  dplyr::arrange(dplyr::desc(cells))

readr::write_csv(qc_pre_summary, file.path(outdir, "qc_pre_summary.csv"))

# 1) counts/filter
msg("# Cell filtering (first pass): nFeature >= ", min_feat, " AND nCount >= ", min_count)
pre_cells <- vapply(seurat_list, ncol, numeric(1))

cf_list <- lapply(
  seurat_list,
  function(obj) {
    keep_cells <- WhichCells(
      obj,
      expression = nFeature_RNA >= min_feat &
        nCount_RNA   >= min_count
    )
    if (length(keep_cells) == 0) {
      msg("  ! Sample dropped by counts/features (0 cells): ", unique(obj$sample_id)[1])
      return(NULL)
    }
    subset(obj, cells = keep_cells)
  }
)

keep_idx <- !vapply(cf_list, is.null, logical(1))
if (!all(keep_idx)) {
  msg("# Samples removed (counts/features -> 0 cells): ",
      paste(names(cf_list)[!keep_idx], collapse = ", "))
}
seurat_list <- cf_list[keep_idx]
post_cf_cells <- vapply(seurat_list, ncol, numeric(1))
msg("# After counts/features filter - samples: ", length(seurat_list),
    " | total cells: ", sum(post_cf_cells))

# 2) scDblFinder

if (isTRUE(run_doublets)) {
  msg("# Running scDblFinder on counts/features-filtered cells...")
  sample_names <- names(seurat_list)
  
  dbl_list <- lapply(seq_along(seurat_list), function(i) {
    obj <- seurat_list[[i]]
    nm  <- sample_names[i]
    
    ok <- TRUE
    sce <- tryCatch(as.SingleCellExperiment(obj), error = function(e) { ok <<- FALSE; NULL })
    if (!ok || is.null(sce)) {
      msg("  ! as.SingleCellExperiment failed: ", nm, " (skipping)")
      return(list(name = nm, obj = obj, dbl_meta = NULL))
    }
    
    ok <- TRUE
    sce <- tryCatch(scDblFinder(sce), error = function(e) { ok <<- FALSE; sce })
    if (!ok) {
      msg("  ! scDblFinder failed: ", nm, " (keeping as-is)")
      return(list(name = nm, obj = obj, dbl_meta = NULL))
    }
    
    dmeta <- as.data.frame(SummarizedExperiment::colData(sce))
    dmeta$barcode <- rownames(dmeta)
    keep_cells <- rownames(obj@meta.data) %in% dmeta$barcode[dmeta$scDblFinder.class == "singlet"]
    obj2 <- subset(obj, cells = rownames(obj@meta.data)[keep_cells])
    list(name = nm, obj = obj2, dbl_meta = dmeta)
  })
  
  seurat_list <- setNames(lapply(dbl_list, `[[`, "obj"), sapply(dbl_list, `[[`, "name"))
  
  
  # Csv cell doublet calls for plotting later
  
  dbl_cells <- dplyr::bind_rows(lapply(dbl_list, function(x) {
    if (is.null(x$dbl_meta)) return(NULL)
    
    d <- x$dbl_meta
    
    # Ensure barcode column
    if (!"barcode" %in% colnames(d)) {
      d$barcode <- rownames(d)
    }
    
    # Keep only useful columns (class is required; others are optional)
    keep_cols <- intersect(
      c("barcode", "scDblFinder.class", "scDblFinder.score", "scDblFinder.weight"),
      colnames(d)
    )
    
    dplyr::as_tibble(d[, keep_cols, drop = FALSE]) %>%
      dplyr::mutate(sample_id = x$name, .before = 1) %>%
      dplyr::rename(cell_barcode = barcode)
  }))
  
  if (!is.null(dbl_cells) && nrow(dbl_cells) > 0) {
    readr::write_csv(dbl_cells, file.path(outdir, "doublet_cells.csv"))
    msg("# Saved: doublet_cells.csv (per-cell scDblFinder calls)")
  } else {
    msg("# NOTE: doublet_cells.csv not written (no per-cell dbl_meta available)")
  }
  
  
  #
  dbl_report <- do.call(rbind, lapply(dbl_list, function(x) {
    if (is.null(x$dbl_meta)) {
      data.frame(sample = x$name, singlets = ncol(x$obj), doublets = NA_integer_, frac_doublet = NA_real_)
    } else {
      tab <- table(x$dbl_meta$scDblFinder.class)
      sing <- if ("singlet" %in% names(tab)) as.integer(tab[["singlet"]]) else 0L
      doub <- if ("doublet" %in% names(tab)) as.integer(tab[["doublet"]]) else 0L
      data.frame(sample = x$name, singlets = sing, doublets = doub,
                 frac_doublet = if ((sing + doub) > 0) doub/(sing + doub) else NA_real_)
    }
  }))
  
  msg("# Doublet summary (head):")
  print(head(dbl_report[order(-dbl_report$frac_doublet), ], 10))
  readr::write_csv(dbl_report, file.path(outdir, "doublet_summary.csv"))
} else {
  msg("# scDblFinder skipped (run_doublets = FALSE)")
}

post_dbl_cells <- vapply(seurat_list, ncol, numeric(1))
msg("# After doublet removal - samples: ", length(seurat_list),
    " | total cells: ", sum(post_dbl_cells))

# 3) mt %

msg("# Mito filter (final step): percent.mt <= ", mt_max)
filt_mt_list <- lapply(
  seurat_list,
  function(obj) {
    keep_cells <- WhichCells(obj, expression = percent.mt <= mt_max)
    if (length(keep_cells) == 0) {
      msg("  ! Sample dropped by mito filter (0 cells): ", unique(obj$sample_id)[1])
      return(NULL)
    }
    subset(obj, cells = keep_cells)
  }
)

keep_idx <- !vapply(filt_mt_list, is.null, logical(1))
if (!all(keep_idx)) {
  msg("# Samples removed by mito filter: ",
      paste(names(filt_mt_list)[!keep_idx], collapse = ", "))
}
seurat_list <- filt_mt_list[keep_idx]

post_mt_cells <- vapply(seurat_list, ncol, numeric(1))
msg("# After mito filter - samples: ", length(seurat_list),
    " | total cells: ", sum(post_mt_cells))

# Post-filter summary 
post_summary <- purrr::map2_dfr(
  seurat_list, names(seurat_list),
  function(o, nm) {
    tibble::tibble(
      sample   = nm,
      cells    = ncol(o),
      mean_mt  = mean(o$percent.mt, na.rm = TRUE),
      mean_nF  = mean(o$nFeature_RNA, na.rm = TRUE),
      mean_nC  = mean(o$nCount_RNA,   na.rm = TRUE)
    )
  }
) %>% dplyr::arrange(dplyr::desc(cells))



msg("# Post-filter summary (head):")
print(head(post_summary, 10))

# outputs
readr::write_csv(post_summary, file.path(outdir, "qc_post_summary.csv"))
saveRDS(seurat_list, file.path(outdir, "seurat_list_filtered.rds"))

summary_qc_cells <- tibble::tibble(
  metric = c(
    "cells_pre",
    "cells_after_counts_features",
    "cells_after_doublets",
    "cells_after_mt"
  ),
  value = c(
    sum(pre_cells, na.rm = TRUE),
    sum(post_cf_cells, na.rm = TRUE),
    sum(post_dbl_cells, na.rm = TRUE),
    sum(post_mt_cells, na.rm = TRUE)
  )
)

readr::write_csv(summary_qc_cells, file.path(outdir, "qc_summary_cells_sequential.csv"))

msg("# Saved: qc_summary_cells_sequential.csv")

msg("# DONE. Saved: qc_post_summary.csv, doublet_summary.csv (if run), seurat_list_filtered.rds")
