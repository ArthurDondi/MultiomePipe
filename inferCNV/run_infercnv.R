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

# in a fresh R session, before anything touches Python:
# library(reticulate)
# py_require("leidenalg")
# py_module_available("leidenalg")   # TRUE — this triggers the env build
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(infercnv)
})

library(reticulate)
reticulate::py_require("leidenalg")

# Avoids an hclust error inferCNV can throw in analysis_mode="subclusters".
options(scipen = 100)

# =============================== USER SETTINGS ================================
# Folder produced by export_for_infercnv.py (copy it next to RStudio if needed).
INPUT_DIR  <- "/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/QC/RNA/Merged/inferCNV/input"
# Where inferCNV writes its results (one sub-folder per sample in per-sample mode).
OUTPUT_DIR <- "/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/QC/RNA/Merged/inferCNV/results"

LOG_RUN <- TRUE               # tee inferCNV's own log lines + this script's summaries to a
# timestamped file in OUTPUT_DIR. It CANNOT capture R warnings, other
# packages' messages (Seurat), or C++ progress bars -- for a COMPLETE
# transcript use OS-level capture instead: an sbatch `#SBATCH --output`,
# or from a shell `Rscript run_infercnv.R 2>&1 | tee full.log`.

MODE <- "per_sample"          # "per_sample" (one run per sample) or "joint" (one run)

# Which samples to process. Empty vector = ALL samples found in metadata.tsv.
# In per_sample mode these are looped one run each; in joint mode they are the
# cells pooled into the single run. e.g. SAMPLES <- c("BMO-SKNBE2c", "Fetahu_M1")
SAMPLES <- c("BMO-STANB8", "BMO-STANB12", "BMO-SKNSH", "BMO-SKNBE2c")

# Cell types to DROP entirely before inferCNV (neither reference nor observation),
# e.g. junk / non-cell populations. Must match metadata.tsv$celltype. Empty = keep all.
EXCLUDE_CELL_TYPES <- c("damaged cells", "doublets", "unknown",
                        "neutrophil", "neutr. prog.", "eo/baso", "mast")

# Cell types treated as the NORMAL / diploid reference. Must match the labels in
# metadata.tsv$celltype (the exporter prints the available labels; adjust here).
NORMAL_CELL_TYPES <- c("HSC/MPP","LMPP","ELP","MEP","GMP", "MK", "PS",
                       "mono. prog.","monocyte","mon. mac.","DC",
                       "plasmacytoid dendritic cell","T cell", "hematopoietic precursor cell",
                       "natural killer cell","B cell","erythroid lineage cell",
                       "endothelial rest", "mesenchymal", "epithelial", "pre HE")


# Merge all reference cell types into ONE "reference" group. Strongly recommended
# with USE_HMM = TRUE: inferCNV's HMM h-spike step crashes ("'x' must be an array
# of at least two dimensions") when there are many small reference groups.
# Observations keep their own per-cell-type labels either way.
COLLAPSE_REFERENCE <- TRUE

# inferCNV run parameters
CUTOFF        <- 0.1          # 0.1 for 10x (droplet); use 1 for Smart-seq2
ANALYSIS_MODE <- "subclusters"    # "samples" (per-group calls) or "subclusters" (per-subclone calls)
# Subcluster controls, used only when ANALYSIS_MODE = "subclusters". inferCNV has
# NO hard "max N subclusters" knob; a LOW leiden_resolution biases toward few
# (aim for <=2, then check .observation_groupings.txt). References are never
# subclustered into the CNV predictions. For a GUARANTEED cap, pre-define groups.
TUMOR_SUBCLUSTER_PARTITION <- "leiden"   # "leiden" or "random_trees"
LEIDEN_FUNCTION            <- "modularity" # "modularity" or "CPM" (inferCNV default). modularity +
# a low resolution + large k_nn tends to give fewer clusters.
LEIDEN_RESOLUTION          <- 0.1        # lower -> fewer subclusters; "auto" = inferCNV default
K_NN                       <- 50         # neighbours in the leiden graph (default 20); larger -> fewer
TUMOR_SUBCLUSTER_PVAL      <- 0.1        # random_trees only: lower -> fewer splits

# How to get subclones:
#   "manual" (recommended) -> FAST 2-pass: pass 1 runs "samples" mode WITHOUT the
#      HMM to get the CNV matrix, cuts each tumour group into <= MANUAL_K subclones
#      by hierarchical clustering ON THAT CNV, then pass 2 runs "samples" mode WITH
#      the HMM on those labels. CNV-based, hard <=MANUAL_K cap, reference never
#      subclustered, and per-subclone `hmm_mode-samples...pred_cnv_regions.dat`.
#   "single" -> one inferCNV run using ANALYSIS_MODE. ANALYSIS_MODE="subclusters"
#      is inferCNV's native leiden, which can take MANY HOURS on large/clonal
#      groups (it also subclusters the reference) -- avoid under a session limit.
SUBCLONE_STRATEGY <- "single"
MAX_SUBCLONES     <- 2        # UPPER BOUND per tumour group; the actual number (1..MAX) is chosen
# from the data, so a monoclonal sample stays a single clone.
SUBCLONE_MIN_SIL  <- 0.20     # accept a split only if its mean silhouette >= this (how real the
# structure is). Lower -> split more readily; raise -> demand cleaner
# subclones. Check the logged silhouette values to calibrate.
SUBCLONE_MIN_CELLS <- 50      # do not split a tumour group smaller than this
USE_HMM       <- TRUE         # HMM CNV-state calling -> binarized predictions (needs JAGS)
HMM_TYPE      <- "i6"         # "i6" = 6 CN states (finer; better for MYCN-amp levels) but builds the
print(HMM_TYPE)
# fragile hidden spike-in. "i3" = 3 states (loss/neutral/gain), defined
# from the reference mean+-sd, and SKIPS the spike -> switch to "i3" if a
# sample dies at "-hspike modeling ... 'x' must be an array". Use one
# type for the whole cohort so samples stay comparable.
DENOISE       <- TRUE
NUM_THREADS   <- 4
SAVE_PRELIM_PLOT <- TRUE      # also write the pre-denoising heatmap (infercnv.preliminary.png)

# Guards / options
MIN_CELLS_PER_RUN <- 20       # skip a sample with fewer cells than this
MIN_REF_CELLS     <- 10       # only when COLLAPSE_REFERENCE = FALSE: min cells for a normal type
# to be its OWN reference group (tiny separate groups break the
# h-spike). Ignored when collapsing, where every present normal
# type is pooled regardless of size.
SUBSAMPLE_MAX     <- 1000     # QUERY cap: max cells per OBSERVATION (tumour) cell type per run,
# e.g. 1000. NULL = keep all tumour cells. Applies per cell type.
REF_SUBSAMPLE     <- 500      # REFERENCE cap: pool ALL reference (normal) cells together and randomly
# keep this many total (not per type) -- a stable diploid baseline needs
# a total, not per-type balance. NULL = keep all reference cells.
SEED              <- 1
# =============================================================================

set.seed(SEED)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
GENE_ORDER <- file.path(INPUT_DIR, "gene_ordering.tsv")
stopifnot(file.exists(GENE_ORDER))

# ---- run log ----------------------------------------------------------------
# inferCNV logs through futile.logger; teeing its ROOT logger sends every
# INFO/WARN/timing line to console AND the log file. loginfo() routes this
# script's own summaries through the same logger, so the log is complete.
if (LOG_RUN) {
  LOG_FILE <- file.path(OUTPUT_DIR, format(Sys.time(), "run_infercnv_%Y%m%d_%H%M%S.log"))
  futile.logger::flog.appender(futile.logger::appender.tee(LOG_FILE))
}
loginfo <- function(txt) if (LOG_RUN) futile.logger::flog.info("%s", txt) else message(txt)
if (LOG_RUN) loginfo(sprintf("Logging run to %s", LOG_FILE))

# Thread oversubscription is a classic cause of an absurdly slow run: if
# NUM_THREADS exceeds the CPUs SLURM allocated, inferCNV's threads (parallelDist,
# Seurat) fight over the same core(s) and crawl. Warn loudly when they mismatch.
.slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (!is.na(.slurm_cpus) && NUM_THREADS > .slurm_cpus)
  loginfo(sprintf(paste0("WARNING: NUM_THREADS = %d but SLURM --cpus-per-task = %d -> ",
                         "thread oversubscription (very slow). Set them equal: bump ",
                         "--cpus-per-task to %d, or set NUM_THREADS <- %d."),
                  NUM_THREADS, .slurm_cpus, NUM_THREADS, .slurm_cpus))

# ---- load the exported cohort once ------------------------------------------
loginfo("Loading counts matrix ...")
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

if (length(EXCLUDE_CELL_TYPES) > 0) {            # drop junk populations up front
  drop <- meta$celltype %in% EXCLUDE_CELL_TYPES
  loginfo(sprintf("Excluding %d cells of type(s): %s",
                  sum(drop), paste(EXCLUDE_CELL_TYPES, collapse = ", ")))
  meta <- meta[!drop, , drop = FALSE]
}
# Restrict to the selected samples up front (per_sample/joint both honour SAMPLES),
# then shrink the RESIDENT matrix to only the cells we will actually use. Holding
# the whole cohort in RAM for the entire loop is the main memory overhead per
# sample; excess resident memory is the likeliest reason a step as light as
# clustering a few thousand cells crawls (the machine swaps).
if (length(SAMPLES) > 0 && "sample" %in% colnames(meta)) {
  bad <- setdiff(SAMPLES, unique(meta$sample))
  if (length(bad) > 0) stop("SAMPLES not found in metadata: ", paste(bad, collapse = ", "))
  meta <- meta[meta$sample %in% SAMPLES, , drop = FALSE]
}
counts <- counts[, meta$cell, drop = FALSE]
invisible(gc())
loginfo(sprintf("Loaded %d genes x %d cells; %d samples (resident after subsetting).",
                nrow(counts), ncol(counts),
                if ("sample" %in% colnames(meta)) length(unique(meta$sample)) else NA))

# ---- helpers ----------------------------------------------------------------
subsample_cells <- function(cells, groups, max_per_group) {
  if (is.null(max_per_group)) return(cells)
  unlist(lapply(split(cells, groups), function(cc)
    if (length(cc) > max_per_group) sample(cc, max_per_group) else cc),
    use.names = FALSE)
}

# Build an inferCNV object for `cells` (grouped by `grp`, baseline = `ref_group`)
# and run one pass; returns the run infercnv_obj (its @expr.data is the CNV matrix).
.infercnv_pass <- function(cells, grp, ref_group, out_dir, hmm, analysis_mode,
                           denoise, no_prelim) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ann_path <- file.path(out_dir, "annotations.tsv")
  write.table(data.frame(cell = cells, group = grp[cells]),
              ann_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts[, cells, drop = FALSE],
    annotations_file  = ann_path,
    gene_order_file   = GENE_ORDER,
    ref_group_names   = if (length(ref_group) > 0) ref_group else NULL)
  # NB: CreateInfercnvObject excludes chrX/chrY/chrM by default (chr_exclude).
  # Drop any arg not in THIS inferCNV version's signature so version differences
  # (e.g. leiden_resolution) can never abort the run.
  run_args <- list(
    obj, cutoff = CUTOFF, out_dir = out_dir, cluster_by_groups = TRUE,
    denoise = denoise, HMM = hmm, HMM_type = HMM_TYPE, analysis_mode = analysis_mode,
    num_threads = NUM_THREADS, resume_mode = TRUE, no_prelim_plot = no_prelim,
    output_format = "png",
    tumor_subcluster_partition_method = TUMOR_SUBCLUSTER_PARTITION,
    tumor_subcluster_pval = TUMOR_SUBCLUSTER_PVAL, leiden_resolution = LEIDEN_RESOLUTION,
    leiden_function = LEIDEN_FUNCTION, k_nn = K_NN)
  keep <- names(run_args) == "" | names(run_args) %in% names(formals(infercnv::run))
  do.call(infercnv::run, run_args[keep])
}

# Mean silhouette width of clustering `cl` over distances `d` (how real the
# structure is; near 0 means one blob). Uses the recommended `cluster` package
# when available, else a base-R fallback.
.mean_sil <- function(d, cl) {
  if (length(unique(cl)) < 2) return(-Inf)
  if (requireNamespace("cluster", quietly = TRUE))
    return(mean(cluster::silhouette(cl, d)[, "sil_width"]))
  D <- as.matrix(d); ks <- unique(cl)
  mean(vapply(seq_along(cl), function(i) {
    same <- which(cl == cl[i] & seq_along(cl) != i)
    a <- if (length(same)) mean(D[i, same]) else 0
    b <- min(vapply(setdiff(ks, cl[i]), function(k) mean(D[i, cl == k]), numeric(1)))
    if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }, numeric(1)))
}

# Pick the number of subclones 1..max_k from CNV distances: the k with the best
# mean silhouette, or 1 if even the best is below min_sil. Returns list(k, cl, sil).
.choose_k <- function(d, max_k, min_sil) {
  n <- attr(d, "Size")
  if (max_k < 2 || n < 4) return(list(k = 1L, cl = rep(1L, n), sil = -Inf))
  hc <- hclust(d, method = "ward.D2")
  best <- list(k = 1L, cl = rep(1L, n), sil = -Inf)
  for (k in 2:min(max_k, n - 1L)) {
    cl <- cutree(hc, k = k); s <- .mean_sil(d, cl)
    if (s > best$sil) best <- list(k = k, cl = cl, sil = s)
  }
  if (best$sil < min_sil) list(k = 1L, cl = rep(1L, n), sil = best$sil) else best
}

# Split each observation (tumour) group into its DATA-DRIVEN number of subclones
# (1..max_k) by clustering the inferred-CNV matrix `cnv` (genes x cells). Reference
# groups (`ref_labels`) are never touched. Returns a per-cell label vector.
.split_observations <- function(cnv, grp, ref_labels, max_k, min_cells, min_sil) {
  sub <- grp
  for (g in setdiff(unique(grp), ref_labels)) {
    cc <- intersect(names(grp)[grp == g], colnames(cnv))
    if (length(cc) < max(min_cells, 4L)) next               # too small to split
    r <- .choose_k(dist(t(cnv[, cc, drop = FALSE])), max_k, min_sil)
    if (r$k > 1L) {
      sub[cc] <- paste0(g, "::", r$cl)
      loginfo(sprintf("  '%s': %d subclones (silhouette %.2f): %s", g, r$k, r$sil,
                      paste(sprintf("%s=%d", names(table(r$cl)), as.integer(table(r$cl))),
                            collapse = ", ")))
    } else {
      loginfo(sprintf("  '%s': 1 clone (best silhouette %.2f < %.2f, no split)",
                      g, r$sil, min_sil))
    }
  }
  sub
}

# Run inferCNV on a set of cells; returns TRUE if a run was launched.
run_one <- function(cells, out_dir, label) {
  grp <- setNames(meta[cells, "celltype"], cells)
  
  # Subsample. REFERENCE (normal) cells are POOLED across all normal types and
  # randomly capped at REF_SUBSAMPLE (a diploid baseline needs a stable total, not
  # per-type balance). OBSERVATION/tumour ("query") cells are capped PER CELL TYPE
  # at SUBSAMPLE_MAX so subclones survive. Either cap can be NULL (keep all).
  is_norm   <- grp %in% NORMAL_CELL_TYPES
  ref_cells <- cells[is_norm]
  obs_cells <- cells[!is_norm]
  if (!is.null(REF_SUBSAMPLE) && length(ref_cells) > REF_SUBSAMPLE)
    ref_cells <- sample(ref_cells, REF_SUBSAMPLE)
  if (!is.null(SUBSAMPLE_MAX))
    obs_cells <- subsample_cells(obs_cells, grp[obs_cells], SUBSAMPLE_MAX)
  cells <- c(ref_cells, obs_cells)
  grp   <- grp[cells]
  
  tab <- sort(table(grp), decreasing = TRUE)
  # reference = normal types present in this run. When kept as SEPARATE groups
  # (COLLAPSE_REFERENCE = FALSE) drop tiny ones (they make inferCNV's h-spike
  # fragile); when collapsed into one group per-type size is irrelevant, so every
  # present normal type contributes — including rare T / NK / B / HSPC.
  ref <- intersect(NORMAL_CELL_TYPES, names(tab))
  if (!COLLAPSE_REFERENCE) ref <- ref[tab[ref] >= MIN_REF_CELLS]
  n_ref <- sum(tab[ref])
  n_obs <- length(cells) - n_ref
  
  loginfo(sprintf("=== %s : %d cells ===", label, length(cells)))
  loginfo(paste(capture.output(print(tab)), collapse = "\n"))
  if (length(ref) > 0) {
    loginfo(sprintf("  reference (normal): %s  [%d cells%s]",
                    paste(ref, collapse = ", "), n_ref,
                    if (COLLAPSE_REFERENCE) " -> merged into one 'reference' group" else ""))
  } else {
    loginfo("  no normal reference cells found -> REFERENCE-LESS run (baseline = average of all cells).")
  }
  if (n_obs == 0) {
    loginfo("  no observation (non-reference) cells -> skipping.")
    return(FALSE)
  }
  
  # Collapse all reference cell types into a single 'reference' group so the HMM
  # h-spike models one large group (observations keep their per-cell-type labels).
  if (COLLAPSE_REFERENCE && length(ref) > 0) {
    grp[grp %in% ref] <- "reference"
    ref <- "reference"
  }
  
  if (SUBCLONE_STRATEGY == "manual") {
    # 2-pass: derive <=MANUAL_K CNV subclones per tumour group, then HMM on them.
    loginfo("  manual subclones - pass 1: CNV matrix (samples mode, no HMM) ...")
    obj1 <- .infercnv_pass(cells, grp, ref, file.path(out_dir, "pass1_cnv"),
                           hmm = FALSE, analysis_mode = "samples",
                           denoise = FALSE, no_prelim = TRUE)
    grp <- .split_observations(obj1@expr.data, grp, ref, MAX_SUBCLONES,
                               SUBCLONE_MIN_CELLS, SUBCLONE_MIN_SIL)
    rm(obj1); gc()
    loginfo("  manual subclones - pass 2: HMM on subclone labels (samples mode) ...")
    .infercnv_pass(cells, grp, ref, out_dir, hmm = USE_HMM, analysis_mode = "samples",
                   denoise = DENOISE, no_prelim = !SAVE_PRELIM_PLOT)
  } else {
    .infercnv_pass(cells, grp, ref, out_dir, hmm = USE_HMM, analysis_mode = ANALYSIS_MODE,
                   denoise = DENOISE, no_prelim = !SAVE_PRELIM_PLOT)
  }
  TRUE
}

# Restrict `all_samples` to the SAMPLES whitelist (empty = keep all); errors
# loudly on an unknown name so a typo can't silently skip a sample.
select_samples <- function(all_samples) {
  if (length(SAMPLES) == 0) return(all_samples)
  missing <- setdiff(SAMPLES, all_samples)
  if (length(missing) > 0)
    stop("SAMPLES not found in metadata: ", paste(missing, collapse = ", "),
         "\n  available: ", paste(all_samples, collapse = ", "))
  intersect(all_samples, SAMPLES)   # keeps all_samples' (sorted) order
}

# ---- dispatch ---------------------------------------------------------------
if (MODE == "per_sample") {
  stopifnot("sample" %in% colnames(meta))
  samples <- select_samples(sort(unique(meta$sample)))
  loginfo(sprintf("Per-sample mode: %d sample(s)%s.", length(samples),
                  if (length(SAMPLES) > 0) " (selected)" else ""))
  for (s in samples) {
    cells <- meta$cell[meta$sample == s]
    if (length(cells) < MIN_CELLS_PER_RUN) {
      loginfo(sprintf("=== %s : %d cells -> below MIN_CELLS_PER_RUN, skipping ===",
                      s, length(cells)))
      next
    }
    out_dir <- file.path(OUTPUT_DIR, gsub("[^A-Za-z0-9._-]", "_", s))
    # Isolate each sample: log a failure and keep going instead of aborting the loop.
    tryCatch(run_one(cells, out_dir, s),
             error = function(e) loginfo(sprintf("!! sample %s FAILED: %s",
                                                 s, conditionMessage(e))))
    invisible(gc())   # release this sample's matrices before the next
  }
} else if (MODE == "joint") {
  cells <- meta$cell
  label <- "joint (all samples)"
  if (length(SAMPLES) > 0) {           # pool only the selected samples
    stopifnot("sample" %in% colnames(meta))
    sel <- select_samples(sort(unique(meta$sample)))
    cells <- meta$cell[meta$sample %in% sel]
    label <- sprintf("joint (%d selected samples)", length(sel))
  }
  run_one(cells, file.path(OUTPUT_DIR, "joint"), label)
} else {
  stop("MODE must be 'per_sample' or 'joint'.")
}

loginfo(paste0("Done. Key outputs per run folder: infercnv.png (heatmap), ",
               "infercnv.observations.txt (per-cell residual expression), and ",
               "HMM_CNV_predictions.* (CNV-state calls) when HMM = TRUE.",
               if (LOG_RUN) sprintf(" Log: %s", LOG_FILE) else ""))