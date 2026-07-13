#!/usr/bin/env Rscript
# =============================================================================
# Run inferCNV on the MultiomePipe batch-corrected cohort (exported by
# inferCNV/export_for_infercnv.py) from R / RStudio.
#
# Neuroblastoma model: the malignant "Neuroblastoma cell" population is compared
# against normal haematopoietic / immune cells (T / NK / B / myeloid / erythroid
# / pDC / HSPC) used as the diploid REFERENCE. inferCNV then infers large-scale
# copy-number gains/losses (e.g. MYCN/2p, 17q gain, 1p / 11q loss) from average
# expression along each chromosome.
#
# With 28 samples the default is a PER-SAMPLE loop (MODE = "per_sample"): every
# sample is run on its own, using that sample's own normal cells as reference.
# This controls for cross-sample / ambient batch effects and keeps each run
# small. Pure cell-line samples that contain no normal cells fall back to a
# reference-less run (baseline = average of all cells in that sample). Set
# MODE = "joint" to run all 28 samples together instead.
#
# ---- one-time install (in your RStudio R) --------------------------------
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("infercnv")        # pulls Matrix, etc.
#   # inferCNV also needs JAGS for the HMM: https://mcmc-jags.sourceforge.io/
#   # (conda: `conda install -c conda-forge jags`; macOS: `brew install jags`)
# --------------------------------------------------------------------------
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(infercnv)
})

# =============================== USER SETTINGS ================================
# Folder produced by export_for_infercnv.py (copy it next to RStudio if needed).
INPUT_DIR  <- "/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/inferCNV/input"
# Where inferCNV writes its results (one sub-folder per sample in per-sample mode).
OUTPUT_DIR <- "/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/inferCNV/results"

MODE <- "per_sample"          # "per_sample" (28 runs) or "joint" (one run)

# Cell types treated as the NORMAL / diploid reference. Must match the labels in
# metadata.tsv$celltype (the exporter prints the available labels; adjust here).
NORMAL_CELL_TYPES <- c(
  "T cell", "NK cell", "B cell", "Myeloid cell",
  "Erythroid lineage cell", "Plasmacytoid dendritic cell",
  "Hematopoietic precursor cell"
)

# inferCNV run parameters
CUTOFF        <- 0.1          # 0.1 for 10x (droplet); use 1 for Smart-seq2
ANALYSIS_MODE <- "samples"    # "samples" (fast) or "subclusters" (tumour subclones; slower, needs leiden)
USE_HMM       <- TRUE         # HMM CNV-state calling (needs JAGS)
DENOISE       <- TRUE
NUM_THREADS   <- 4

# Guards / options
MIN_CELLS_PER_RUN <- 20       # skip a sample with fewer cells than this
MIN_REF_CELLS     <- 10       # a normal type needs >= this many cells to be used as reference
SUBSAMPLE_MAX     <- NULL     # e.g. 2000 to cap cells per group per run for speed; NULL = no subsampling
SEED              <- 1
# =============================================================================

set.seed(SEED)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
GENE_ORDER <- file.path(INPUT_DIR, "gene_ordering.tsv")
stopifnot(file.exists(GENE_ORDER))

# ---- load the exported cohort once ------------------------------------------
message("Loading counts matrix ...")
counts <- readMM(gzfile(file.path(INPUT_DIR, "counts.mtx.gz")))   # genes x cells
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- readLines(file.path(INPUT_DIR, "genes.tsv"))
colnames(counts) <- readLines(file.path(INPUT_DIR, "barcodes.tsv"))

meta <- read.delim(file.path(INPUT_DIR, "metadata.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE,
                   quote = "", comment.char = "")
stopifnot(all(c("cell", "celltype") %in% colnames(meta)))
rownames(meta) <- meta$cell
meta <- meta[colnames(counts), , drop = FALSE]   # align to matrix column order
message(sprintf("Loaded %d genes x %d cells; %d samples.",
                nrow(counts), ncol(counts),
                if ("sample" %in% colnames(meta)) length(unique(meta$sample)) else NA))

# ---- helpers ----------------------------------------------------------------
subsample_cells <- function(cells, groups, max_per_group) {
  if (is.null(max_per_group)) return(cells)
  unlist(lapply(split(cells, groups), function(cc)
    if (length(cc) > max_per_group) sample(cc, max_per_group) else cc),
    use.names = FALSE)
}

# Run inferCNV on a set of cells; returns TRUE if a run was launched.
run_one <- function(cells, out_dir, label) {
  grp <- setNames(meta[cells, "celltype"], cells)

  if (!is.null(SUBSAMPLE_MAX)) {
    cells <- subsample_cells(cells, grp[cells], SUBSAMPLE_MAX)
    grp   <- grp[cells]
  }

  tab <- sort(table(grp), decreasing = TRUE)
  # reference = normal types that are present AND have enough cells
  ref <- intersect(NORMAL_CELL_TYPES, names(tab)[tab >= MIN_REF_CELLS])
  n_ref <- sum(tab[ref])
  n_obs <- length(cells) - n_ref

  message(sprintf("\n=== %s : %d cells ===", label, length(cells)))
  print(tab)
  if (length(ref) > 0) {
    message(sprintf("  reference (normal) groups: %s  [%d cells]",
                    paste(ref, collapse = ", "), n_ref))
  } else {
    message("  no normal reference cells found -> REFERENCE-LESS run ",
            "(baseline = average of all cells).")
  }
  if (n_obs == 0) {
    message("  no observation (non-reference) cells -> skipping.")
    return(FALSE)
  }

  # inferCNV annotations file for this run: <cell> <TAB> <group>, no header.
  ann_path <- file.path(out_dir, "annotations.tsv")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(data.frame(cell = cells, group = grp[cells]),
              ann_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts[, cells, drop = FALSE],
    annotations_file  = ann_path,
    gene_order_file   = GENE_ORDER,
    ref_group_names   = if (length(ref) > 0) ref else NULL
    # NB: CreateInfercnvObject excludes chrX/chrY/chrM by default (chr_exclude);
    # pass chr_exclude = c("chrM") here if you want to keep the sex chromosomes.
  )

  infercnv::run(
    obj,
    cutoff            = CUTOFF,
    out_dir           = out_dir,
    cluster_by_groups = TRUE,          # cluster observations by their cell-type group
    denoise           = DENOISE,
    HMM               = USE_HMM,
    analysis_mode     = ANALYSIS_MODE,
    num_threads       = NUM_THREADS,
    resume_mode       = TRUE,          # re-runs pick up where they left off
    no_prelim_plot    = TRUE,
    output_format     = "png"          # heatmap format; "pdf" for a vector figure
  )
  TRUE
}

# ---- dispatch ---------------------------------------------------------------
if (MODE == "per_sample") {
  stopifnot("sample" %in% colnames(meta))
  samples <- sort(unique(meta$sample))
  message(sprintf("Per-sample mode: %d samples.", length(samples)))
  for (s in samples) {
    cells <- meta$cell[meta$sample == s]
    if (length(cells) < MIN_CELLS_PER_RUN) {
      message(sprintf("\n=== %s : %d cells -> below MIN_CELLS_PER_RUN, skipping ===",
                      s, length(cells)))
      next
    }
    out_dir <- file.path(OUTPUT_DIR, gsub("[^A-Za-z0-9._-]", "_", s))
    run_one(cells, out_dir, s)
  }
} else if (MODE == "joint") {
  run_one(meta$cell, file.path(OUTPUT_DIR, "joint"), "joint (all samples)")
} else {
  stop("MODE must be 'per_sample' or 'joint'.")
}

message("\nDone. Key outputs per run folder: infercnv.png (heatmap), ",
        "infercnv.observations.txt (per-cell residual expression), and ",
        "HMM_CNV_predictions.* (CNV-state calls) when HMM = TRUE.")
