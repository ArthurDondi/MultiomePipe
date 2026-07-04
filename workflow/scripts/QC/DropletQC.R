#!/usr/bin/env Rscript
# DropletQC.R
# 2-step DropletQC: identify_empty_drops() then identify_damaged_cells()
# on ALL raw barcodes. Supports scRNA and snRNA assay types.
# Outputs a comparison table against CellRanger filtered barcodes.

suppressPackageStartupMessages({
    library(argparse)
    library(Seurat)
    library(Matrix)
    # DropletQC is installed at conda-env build time by
    # workflow/envs/dropletqc.post-deploy.sh, so it (and its parallel future
    # workers) can attach the namespace without a fragile runtime install.
    library(DropletQC)
})

parser <- ArgumentParser(description = "DropletQC 2-step damaged-cell detection")
parser$add_argument("--raw_h5",      required = TRUE,
                    help = "Cell Ranger raw_feature_bc_matrix.h5 (all barcodes)")
parser$add_argument("--filtered_h5", required = TRUE,
                    help = "Cell Ranger filtered_feature_bc_matrix.h5 (for comparison)")
parser$add_argument("--bam",         required = TRUE,
                    help = "BAM file (possorted_genome_bam.bam)")
parser$add_argument("--output_dir",  required = TRUE,
                    help = "Directory for output files")
parser$add_argument("--assay_type",  default = "scrna",
                    help = "Assay type: 'scrna' or 'snrna' (default: scrna)")
parser$add_argument("--seurat_resolution", type = "double", default = 0.5,
                    help = "Seurat clustering resolution for cell-type labels (default: 0.5)")
parser$add_argument("--min_umi_for_nf", type = "integer", default = 100,
                    help = "Min UMI count to attempt nuclear fraction computation (default: 100)")
parser$add_argument("--cores", type = "integer", default = 1,
                    help = paste("Worker processes for nuclear_fraction_tags()",
                                 "(default: 1). Should match the rule's allocated",
                                 "threads to avoid oversubscribing the node."))
args <- parser$parse_args()

if (!args$assay_type %in% c("scrna", "snrna")) {
    stop("--assay_type must be 'scrna' or 'snrna', got: ", args$assay_type)
}
if (!file.exists(args$bam)) {
    stop("BAM file not found: ", args$bam)
}

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

load_gex <- function(h5_path) {
    mat <- Read10X_h5(h5_path, use.names = TRUE)
    if (is.list(mat)) mat <- mat[["Gene Expression"]]
    mat
}

# ── 1. Load matrices ──────────────────────────────────────────────────────────

message("[1/6] Loading raw matrix: ", args$raw_h5)
raw_mat <- load_gex(args$raw_h5)
message("  Raw barcodes: ", ncol(raw_mat))

message("[1/6] Loading filtered matrix: ", args$filtered_h5)
filt_mat <- load_gex(args$filtered_h5)
cellranger_bcs <- colnames(filt_mat)
message("  CellRanger filtered barcodes: ", length(cellranger_bcs))

# ── 2. Nuclear fraction on raw barcodes with sufficient UMI ──────────────────

raw_umi <- Matrix::colSums(raw_mat)

bcs_for_nf <- names(raw_umi[raw_umi >= args$min_umi_for_nf])
message(sprintf("[2/6] Computing nuclear fraction for %d barcodes (UMI >= %d) ...",
                length(bcs_for_nf), args$min_umi_for_nf))

nf_result <- nuclear_fraction_tags(
    bam      = args$bam,
    barcodes = bcs_for_nf,
    cores    = args$cores,
    verbose  = FALSE
)

nf_df <- as.data.frame(nf_result)
colnames(nf_df)[colnames(nf_df) == "nuclear_fraction"] <- "nf_umi"
nf_df$barcode <- rownames(nf_df)

qc <- data.frame(
    barcode    = bcs_for_nf,
    umi_counts = as.integer(raw_umi[bcs_for_nf]),
    stringsAsFactors = FALSE
)
qc <- merge(qc, nf_df[, c("barcode", "nf_umi")], by = "barcode", all.x = TRUE)
rownames(qc) <- qc$barcode

message(sprintf("  Nuclear fraction computed for %d barcodes; %d NAs",
                nrow(qc), sum(is.na(qc$nf_umi))))

# ── 3. Step 1 - identify_empty_drops ─────────────────────────────────────────

message("[3/6] Running identify_empty_drops() ...")
qc_complete <- qc[!is.na(qc$nf_umi), ]

# Only pass CellRanger filtered barcodes to identify_empty_drops.
# The EM separates cells (high nf) from empty droplets that slipped through
# the CellRanger filter (low nf). Non-CellRanger barcodes are assigned
# empty_droplet by default.
qc_cr <- qc_complete[qc_complete$barcode %in% cellranger_bcs, ]
message(sprintf("  CellRanger barcodes with valid nf: %d / %d",
                nrow(qc_cr), length(cellranger_bcs)))

ed_input <- data.frame(
    nf_umi     = qc_cr$nf_umi,
    umi_counts = qc_cr$umi_counts,
    stringsAsFactors = FALSE
)
ed_result <- identify_empty_drops(
    ed_input,
    nf_rescue    = 0.05,
    umi_rescue   = 1000,
    include_plot = TRUE,
    plot_name    = "identify_empty_drops.png",
    plot_path    = args$output_dir,
    pdf_png      = "png"
)
message("  identify_empty_drops output columns: ", paste(colnames(ed_result), collapse = ", "))
qc_cr$droplet_class <- ed_result[[ncol(ed_result)]]

n_cells_ed <- sum(qc_cr$droplet_class == "cell", na.rm = TRUE)
n_empty_ed <- sum(qc_cr$droplet_class == "empty_droplet", na.rm = TRUE)
message(sprintf("  identify_empty_drops: %d cells, %d empty droplets",
                n_cells_ed, n_empty_ed))

# Propagate droplet_class back to qc_complete; non-CellRanger barcodes = empty_droplet
qc_complete$droplet_class <- "empty_droplet"
rownames(qc_cr) <- qc_cr$barcode
qc_complete[qc_complete$barcode %in% qc_cr$barcode, "droplet_class"] <-
    qc_cr[qc_complete$barcode[qc_complete$barcode %in% qc_cr$barcode], "droplet_class"]

# Fallback: if EM returned degenerate result, use all CellRanger barcodes as cells
if (n_cells_ed == 0 || n_empty_ed == 0) {
    message("  WARNING: identify_empty_drops returned degenerate result (cells=",
            n_cells_ed, ", empty=", n_empty_ed,
            ") - falling back to CellRanger filtered barcodes")
    qc_complete$droplet_class <- ifelse(
        qc_complete$barcode %in% cellranger_bcs, "cell", "empty_droplet"
    )
    n_cells_ed <- sum(qc_complete$droplet_class == "cell")
    message(sprintf("  Fallback: %d cells from CellRanger filtered list", n_cells_ed))
}

# ── 4. Quick Seurat clustering on putative cells for identify_damaged_cells ──

message("[4/6] Quick Seurat clustering for cell-type labels ...")
cell_bcs <- qc_complete$barcode[qc_complete$droplet_class == "cell"]

cell_mat <- raw_mat[, cell_bcs, drop = FALSE]

so <- CreateSeuratObject(counts = cell_mat, min.cells = 0, min.features = 0)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
so <- FindClusters(so, resolution = args$seurat_resolution, verbose = FALSE)

cluster_labels <- paste0("cluster_", as.character(Idents(so)))
names(cluster_labels) <- colnames(so)
message(sprintf("  Found %d clusters", length(unique(cluster_labels))))

# ── 5. Step 2 - identify_damaged_cells ───────────────────────────────────────

message("[5/6] Running identify_damaged_cells() ...")
cell_qc_df <- qc_complete[qc_complete$droplet_class == "cell", ]
cell_qc_df$cell_type <- cluster_labels[cell_qc_df$barcode]
cell_qc_df <- cell_qc_df[!is.na(cell_qc_df$cell_type), ]

dc_input <- data.frame(
    nf_umi        = cell_qc_df$nf_umi,
    umi_counts    = cell_qc_df$umi_counts,
    empty_droplet = cell_qc_df$droplet_class,
    cell_type     = cell_qc_df$cell_type,
    stringsAsFactors = FALSE
)
dc_result <- identify_damaged_cells(dc_input)

# identify_damaged_cells returns list(df = <annotated data frame>, plots = ...)
dc_df <- if (is.data.frame(dc_result)) dc_result else dc_result[["df"]]
if (is.null(dc_df) || !is.data.frame(dc_df)) {
    stop("identify_damaged_cells: could not extract data frame from result (names: ",
         paste(names(dc_result), collapse = ", "), ")")
}
cell_qc_df$damaged_class <- dc_df[["cell_status"]]

n_damaged <- sum(cell_qc_df$damaged_class == "damaged_cell", na.rm = TRUE)
n_intact  <- sum(cell_qc_df$damaged_class == "cell",         na.rm = TRUE)
message(sprintf("  identify_damaged_cells: %d intact cells, %d damaged cells",
                n_intact, n_damaged))

if (args$assay_type == "scrna") {
    message("  [scRNA] Damaged = high nf (cytoplasm leaked, only nucleus retained)")
} else {
    message("  [snRNA] Damaged = low nf among cells (nucleus lost / ambient RNA)")
}

# ── 5b. Combined plot: empty_drops (top) + damaged_cells (bottom) ─────────────

tryCatch({
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(patchwork))

    colours <- c(damaged_cell = "#F8766D", cell = "#00BA38", empty_droplet = "#619CFF")

    # ── top panel: all 3 categories on the full nf × UMI scatter ──────────────
    plot_df <- qc_complete[, c("barcode", "nf_umi", "umi_counts", "droplet_class")]
    plot_df <- merge(plot_df,
                     cell_qc_df[, c("barcode", "damaged_class")],
                     by = "barcode", all.x = TRUE)
    plot_df$category <- plot_df$droplet_class
    plot_df$category[!is.na(plot_df$damaged_class) &
                     plot_df$damaged_class == "damaged_cell"] <- "damaged_cell"
    plot_df$category <- factor(plot_df$category,
                               levels = c("damaged_cell", "cell", "empty_droplet"))

    p_top <- ggplot(plot_df, aes(x = nf_umi, y = log10(umi_counts),
                                 colour = category)) +
        geom_point(size = 0.4, alpha = 0.6) +
        scale_colour_manual(values = colours,
                            breaks = c("damaged_cell", "cell", "empty_droplet")) +
        labs(x = "Nuclear fraction", y = "log10(UMI count)",
             colour = NULL,
             title = "DropletQC: all categories") +
        theme_bw(base_size = 11) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

    # ── bottom panel: cells only, damaged highlighted ──────────────────────────
    cells_df <- cell_qc_df[, c("nf_umi", "umi_counts", "damaged_class")]
    cells_df$damaged_class <- factor(cells_df$damaged_class,
                                     levels = c("damaged_cell", "cell"))

    p_bottom <- ggplot(cells_df, aes(x = nf_umi, y = log10(umi_counts),
                                     colour = damaged_class)) +
        geom_point(size = 0.4, alpha = 0.6) +
        scale_colour_manual(values = c(damaged_cell = "#F8766D", cell = "#00BA38")) +
        labs(x = "Nuclear fraction", y = "log10(UMI count)", colour = NULL,
             title = "identify_damaged_cells: cells only") +
        theme_bw(base_size = 11) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

    combined_path <- file.path(args$output_dir, "dropletqc_combined.png")
    ggsave(combined_path, p_top / p_bottom, width = 6, height = 6, dpi = 150)
    message("  Combined plot written to: ", combined_path)
}, error = function(e) {
    message("  WARNING: could not write combined plot: ", conditionMessage(e))
})

# ── 5c. Cluster plot: nf × UMI coloured by Seurat cluster ─────────────────────
# Damaged cells (high nf, low UMI) tend to collapse into their OWN cluster(s).
# identify_damaged_cells() only flags damage RELATIVE to other cells of the same
# cluster (a within-cluster bimodal split), so a cluster that is uniformly
# damaged is left entirely unflagged. Colouring the nf × UMI scatter by Seurat
# cluster exposes this: if the high-nf / low-UMI population is its own cluster,
# DropletQC's damaged-cell step cannot remove it and a different gate is needed.
tryCatch({
    suppressPackageStartupMessages(library(ggplot2))

    clust_df <- cell_qc_df[, c("nf_umi", "umi_counts", "cell_type")]
    # Order cluster levels numerically (cluster_2 before cluster_10) so the
    # legend and colour assignment are readable and stable across runs.
    clust_levels <- unique(clust_df$cell_type)
    clust_levels <- clust_levels[order(
        suppressWarnings(as.integer(sub("^cluster_", "", clust_levels))),
        clust_levels
    )]
    clust_df$cell_type <- factor(clust_df$cell_type, levels = clust_levels)

    p_clusters <- ggplot(clust_df, aes(x = nf_umi, y = log10(umi_counts),
                                       colour = cell_type)) +
        geom_point(size = 0.4, alpha = 0.6) +
        labs(x = "Nuclear fraction", y = "log10(UMI count)",
             colour = "Seurat cluster",
             title = "DropletQC cells: coloured by Seurat cluster") +
        theme_bw(base_size = 11) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

    clusters_path <- file.path(args$output_dir, "dropletqc_clusters.png")
    ggsave(clusters_path, p_clusters, width = 6.5, height = 4, dpi = 150)
    message("  Cluster plot written to: ", clusters_path)
}, error = function(e) {
    message("  WARNING: could not write cluster plot: ", conditionMessage(e))
})

# ── 6. Build final output & compare with CellRanger ─────────────────────────

message("[6/6] Building outputs and comparing with CellRanger ...")

final <- data.frame(
    barcode    = qc$barcode,
    umi_counts = qc$umi_counts,
    nf_umi     = qc$nf_umi,
    stringsAsFactors = FALSE
)

final <- merge(final,
               qc_complete[, c("barcode", "droplet_class")],
               by = "barcode", all.x = TRUE)

final <- merge(final,
               cell_qc_df[, c("barcode", "damaged_class", "cell_type")],
               by = "barcode", all.x = TRUE)

final$cell_status <- "below_umi_threshold"
final$cell_status[final$droplet_class == "empty_droplet"] <- "empty_droplet"
final$cell_status[final$droplet_class == "cell" & final$damaged_class == "damaged_cell"] <- "damaged_cell"
final$cell_status[final$droplet_class == "cell" & final$damaged_class == "cell"]         <- "cell"

final$in_cellranger_filtered <- final$barcode %in% cellranger_bcs

dropletqc_cells <- final$barcode[final$cell_status == "cell"]

in_both       <- intersect(dropletqc_cells, cellranger_bcs)
only_dqc      <- setdiff(dropletqc_cells, cellranger_bcs)
only_cellrngr <- setdiff(cellranger_bcs,  dropletqc_cells)

message(sprintf("  DropletQC cells:        %d", length(dropletqc_cells)))
message(sprintf("  CellRanger cells:       %d", length(cellranger_bcs)))
message(sprintf("  In both:                %d", length(in_both)))
message(sprintf("  Only in DropletQC:      %d", length(only_dqc)))
message(sprintf("  Only in CellRanger:     %d", length(only_cellrngr)))

# ── Write outputs ─────────────────────────────────────────────────────────────

cell_qc_path <- file.path(args$output_dir, "cell_qc.tsv")
write.table(
    final[, c("barcode", "cell_status", "nf_umi", "umi_counts",
              "droplet_class", "damaged_class", "cell_type",
              "in_cellranger_filtered")],
    cell_qc_path,
    sep = "\t", quote = FALSE, row.names = FALSE
)
message("Written full cell QC to: ", cell_qc_path)

comparison <- data.frame(
    category = c("dropletqc_cells", "cellranger_cells",
                 "in_both", "only_in_dropletqc", "only_in_cellranger"),
    n        = c(length(dropletqc_cells), length(cellranger_bcs),
                 length(in_both), length(only_dqc), length(only_cellrngr))
)
comparison_path <- file.path(args$output_dir, "cellranger_comparison.tsv")
write.table(comparison, comparison_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Written CellRanger comparison to: ", comparison_path)

bc_comparison <- data.frame(
    barcode       = unique(c(dropletqc_cells, cellranger_bcs)),
    in_dropletqc  = unique(c(dropletqc_cells, cellranger_bcs)) %in% dropletqc_cells,
    in_cellranger = unique(c(dropletqc_cells, cellranger_bcs)) %in% cellranger_bcs
)
bc_comparison$agreement <- bc_comparison$in_dropletqc == bc_comparison$in_cellranger
bc_comparison_path <- file.path(args$output_dir, "barcode_comparison.tsv")
write.table(bc_comparison, bc_comparison_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Written per-barcode comparison to: ", bc_comparison_path)

message("DropletQC complete.")