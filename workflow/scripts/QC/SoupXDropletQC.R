#!/usr/bin/env Rscript
# SoupXDropletQC.R
# Ambient RNA correction (SoupX) with optional damaged-cell detection (DropletQC).
#
# Outputs written to --output_dir:
#   corrected/         – 10x-format sparse matrix (SoupX-corrected GEX counts)
#   cell_qc.tsv        – per-cell QC: barcode, soupx_rho [, nf_umi, cell_status]
#   soup_profile.tsv   – ambient gene profile ranked by estimated fraction

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(SoupX)
    library(DropletUtils)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser <- ArgumentParser(description = "SoupX ambient correction + optional DropletQC")
parser$add_argument("--raw_h5",       required = TRUE,
                    help = "Cell Ranger raw_feature_bc_matrix.h5")
parser$add_argument("--filtered_h5",  required = TRUE,
                    help = "Cell Ranger filtered_feature_bc_matrix.h5")
parser$add_argument("--output_dir",   required = TRUE,
                    help = "Directory for all outputs")
parser$add_argument("--bam",          default = NULL,
                    help = "BAM file for DropletQC NF computation (optional)")
parser$add_argument("--contaminant_genes", nargs = "*", default = NULL,
                    help = "Known ambient genes for SoupX (optional)")
parser$add_argument("--resolution",   type = "double", default = 1.0,
                    help = "Seurat clustering resolution for SoupX (default 1.0)")
parser$add_argument("--min_nf_umi",   type = "double", default = 0.6,
                    help = "DropletQC nuclear fraction threshold (default 0.6)")
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(args$output_dir, "corrected"), showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper: load 10x h5, handling Multiome (list) vs single-modal
# ---------------------------------------------------------------------------
load_gex <- function(h5_path) {
    mat <- Read10X_h5(h5_path, use.names = TRUE)
    if (is.list(mat)) {
        message("Multiome h5 detected – extracting 'Gene Expression' slot")
        mat <- mat[["Gene Expression"]]
    }
    mat
}

# ---------------------------------------------------------------------------
# 1. Load raw and filtered matrices
# ---------------------------------------------------------------------------
message("Loading raw matrix: ", args$raw_h5)
raw_mat <- load_gex(args$raw_h5)

message("Loading filtered matrix: ", args$filtered_h5)
filt_mat <- load_gex(args$filtered_h5)

n_cells_filtered <- ncol(filt_mat)
message("Cell Ranger filtered cell count: ", n_cells_filtered)

# ---------------------------------------------------------------------------
# 2. Optional DropletQC – detect damaged cells via nuclear fraction
# ---------------------------------------------------------------------------
cell_qc <- data.frame(barcode = colnames(filt_mat), stringsAsFactors = FALSE)
rownames(cell_qc) <- cell_qc$barcode

run_dropletqc <- !is.null(args$bam) && nchar(args$bam) > 0 && file.exists(args$bam)

if (run_dropletqc) {
    message("Running DropletQC nuclear fraction computation from BAM: ", args$bam)
    suppressPackageStartupMessages(library(DropletQC))
    nf <- nuclear_fraction_tags(
        bam = args$bam,
        barcodes = colnames(filt_mat),
        verbose = FALSE
    )
    nf_df <- as.data.frame(nf)
    nf_df$barcode <- rownames(nf_df)
    cell_qc <- merge(cell_qc, nf_df[, c("barcode", "nf_umi")], by = "barcode", all.x = TRUE)
    rownames(cell_qc) <- cell_qc$barcode

    # Flag damaged cells
    cell_qc$cell_status <- ifelse(
        !is.na(cell_qc$nf_umi) & cell_qc$nf_umi < args$min_nf_umi,
        "damaged_cell",
        "cell"
    )
    n_damaged <- sum(cell_qc$cell_status == "damaged_cell", na.rm = TRUE)
    message(sprintf("DropletQC: %d / %d cells flagged as damaged (nf_umi < %.2f)",
                    n_damaged, nrow(cell_qc), args$min_nf_umi))

    # Remove damaged cells from filtered matrix before SoupX
    keep <- cell_qc$barcode[cell_qc$cell_status == "cell"]
    filt_mat <- filt_mat[, keep, drop = FALSE]
    message("Cells remaining after DropletQC filter: ", ncol(filt_mat))
} else {
    if (!is.null(args$bam) && nchar(args$bam) > 0) {
        warning("BAM file specified but not found: ", args$bam, " – skipping DropletQC")
    }
    cell_qc$cell_status <- "cell"
}

# ---------------------------------------------------------------------------
# 3. Quick Seurat clustering to provide SoupX with cell-type labels
# ---------------------------------------------------------------------------
message("Running quick Seurat clustering for SoupX...")
suppressPackageStartupMessages(library(Seurat))

so <- CreateSeuratObject(counts = filt_mat, min.cells = 3, min.features = 200)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = args$resolution, verbose = FALSE)

clusters <- setNames(as.character(so$seurat_clusters), colnames(so))
message(sprintf("Found %d clusters (resolution = %.1f)", length(unique(clusters)), args$resolution))

# ---------------------------------------------------------------------------
# 4. SoupX ambient RNA correction
# ---------------------------------------------------------------------------
message("Creating SoupX SoupChannel object...")
sc_obj <- SoupChannel(tod = raw_mat, toc = filt_mat)
sc_obj <- setClusters(sc_obj, clusters)

if (!is.null(args$contaminant_genes) && length(args$contaminant_genes) > 0) {
    present_genes <- intersect(args$contaminant_genes, rownames(filt_mat))
    if (length(present_genes) > 0) {
        message("Using user-supplied contaminant genes: ", paste(present_genes, collapse = ", "))
        sc_obj <- setContaminationFraction(sc_obj,
                                           useToEst = present_genes)
    } else {
        warning("None of the provided contaminant genes found in matrix – using autoEstCont")
    }
}

message("Estimating contamination with autoEstCont...")
sc_obj <- tryCatch(
    autoEstCont(sc_obj, verbose = FALSE),
    error = function(e) {
        message("autoEstCont failed (", conditionMessage(e), "), falling back to fixed rho = 0.05")
        setContaminationFraction(sc_obj, 0.05)
    }
)

message("Adjusting counts...")
corrected_mat <- adjustCounts(sc_obj, roundToInt = TRUE, verbose = FALSE)
message(sprintf("SoupX corrected %d genes x %d cells", nrow(corrected_mat), ncol(corrected_mat)))

# ---------------------------------------------------------------------------
# 5. Collect per-cell contamination fraction (rho)
# ---------------------------------------------------------------------------
rho_df <- sc_obj$metaData[, "rho", drop = FALSE]
rho_df$barcode <- rownames(rho_df)
colnames(rho_df)[colnames(rho_df) == "rho"] <- "soupx_rho"
cell_qc <- merge(cell_qc, rho_df, by = "barcode", all.x = TRUE)
rownames(cell_qc) <- cell_qc$barcode

# For cells removed by DropletQC, rho will be NA – mark them explicitly
if ("cell_status" %in% colnames(cell_qc)) {
    cell_qc$soupx_rho[cell_qc$cell_status == "damaged_cell"] <- NA_real_
}

# ---------------------------------------------------------------------------
# 6. Write corrected matrix (10x sparse format)
# ---------------------------------------------------------------------------
message("Writing corrected matrix to: ", file.path(args$output_dir, "corrected"))
write10xCounts(
    path      = file.path(args$output_dir, "corrected"),
    x         = corrected_mat,
    type      = "sparse",
    overwrite = TRUE
)

# ---------------------------------------------------------------------------
# 7. Write per-cell QC table
# ---------------------------------------------------------------------------
cell_qc_out <- cell_qc[, c("barcode", "soupx_rho", "cell_status"), drop = FALSE]
if (run_dropletqc && "nf_umi" %in% colnames(cell_qc)) {
    cell_qc_out <- cbind(cell_qc_out,
                         nf_umi = cell_qc[cell_qc_out$barcode, "nf_umi"])
}
cell_qc_path <- file.path(args$output_dir, "cell_qc.tsv")
write.table(cell_qc_out, cell_qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written cell QC to: ", cell_qc_path)

# ---------------------------------------------------------------------------
# 8. Write ambient gene profile
# ---------------------------------------------------------------------------
soup_prof <- sc_obj$soupProfile[order(-sc_obj$soupProfile$est), ]
soup_prof$gene <- rownames(soup_prof)
soup_path <- file.path(args$output_dir, "soup_profile.tsv")
write.table(soup_prof[, c("gene", "est", "cnts")], soup_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Written soup profile to: ", soup_path)

message("Done.")
