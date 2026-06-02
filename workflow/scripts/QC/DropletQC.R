#!/usr/bin/env Rscript
# DropletQC.R
# Compute nuclear fraction from a BAM file and flag damaged cells.

suppressPackageStartupMessages({
    library(argparse)
    library(DropletQC)
    library(DropletUtils)
})

parser <- ArgumentParser(description = "DropletQC damaged-cell detection")
parser$add_argument("--filtered_h5", required = TRUE,
                    help = "Cell Ranger filtered_feature_bc_matrix.h5")
parser$add_argument("--output_dir", required = TRUE,
                    help = "Directory for DropletQC outputs")
parser$add_argument("--bam", required = TRUE,
                    help = "BAM file for nuclear fraction computation")
parser$add_argument("--min_nf_umi", type = "double", default = 0.6,
                    help = "DropletQC nuclear fraction threshold (default 0.6)")
args <- parser$parse_args()

if (!file.exists(args$bam)) {
    stop("DropletQC BAM file not found: ", args$bam)
}

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

load_gex <- function(h5_path) {
    mat <- Read10X_h5(h5_path, use.names = TRUE)
    if (is.list(mat)) {
        message("Multiome h5 detected – extracting 'Gene Expression' slot")
        mat <- mat[["Gene Expression"]]
    }
    mat
}

message("Loading filtered matrix: ", args$filtered_h5)
filt_mat <- load_gex(args$filtered_h5)

message("Running DropletQC nuclear fraction computation from BAM: ", args$bam)
nf <- nuclear_fraction_tags(
    bam = args$bam,
    barcodes = colnames(filt_mat),
    verbose = FALSE
)

cell_qc <- data.frame(barcode = colnames(filt_mat), stringsAsFactors = FALSE)
rownames(cell_qc) <- cell_qc$barcode

nf_df <- as.data.frame(nf)
nf_df$barcode <- rownames(nf_df)
cell_qc <- merge(cell_qc, nf_df[, c("barcode", "nf_umi")], by = "barcode", all.x = TRUE)
rownames(cell_qc) <- cell_qc$barcode

cell_qc$cell_status <- ifelse(
    !is.na(cell_qc$nf_umi) & cell_qc$nf_umi < args$min_nf_umi,
    "damaged_cell",
    "cell"
)

n_damaged <- sum(cell_qc$cell_status == "damaged_cell", na.rm = TRUE)
message(sprintf("DropletQC: %d / %d cells flagged as damaged (nf_umi < %.2f)",
                n_damaged, nrow(cell_qc), args$min_nf_umi))

cell_qc_path <- file.path(args$output_dir, "cell_qc.tsv")
write.table(
    cell_qc[, c("barcode", "cell_status", "nf_umi"), drop = FALSE],
    cell_qc_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
message("Written cell QC to: ", cell_qc_path)
