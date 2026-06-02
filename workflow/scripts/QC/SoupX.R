#!/usr/bin/env Rscript
# SoupX.R
# Ambient RNA correction (SoupX) after DropletQC damaged-cell filtering.
#
# Outputs written to --output_dir:
#   corrected/         – 10x-format sparse matrix (SoupX-corrected GEX counts)
#   cell_qc.tsv        – per-cell QC: barcode, soupx_rho, cell_status, nf_umi
#   soup_profile.tsv   – ambient gene profile ranked by estimated fraction

suppressPackageStartupMessages({
    library(argparse)
    library(DropletUtils)
    library(SoupX)
})

parser <- ArgumentParser(description = "SoupX ambient correction after DropletQC")
parser$add_argument("--raw_h5", required = TRUE,
                    help = "Cell Ranger raw_feature_bc_matrix.h5")
parser$add_argument("--filtered_h5", required = TRUE,
                    help = "Cell Ranger filtered_feature_bc_matrix.h5")
parser$add_argument("--cell_qc", required = TRUE,
                    help = "DropletQC per-cell QC table")
parser$add_argument("--output_dir", required = TRUE,
                    help = "Directory for all outputs")
parser$add_argument("--contaminant_genes", nargs = "*", default = NULL,
                    help = "Known ambient genes for SoupX (optional)")
parser$add_argument("--resolution", type = "double", default = 1.0,
                    help = "Seurat clustering resolution for SoupX (default 1.0)")
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(args$output_dir, "corrected"), showWarnings = FALSE)

load_gex <- function(h5_path) {
    mat <- Read10X_h5(h5_path, use.names = TRUE)
    if (is.list(mat)) {
        message("Multiome h5 detected – extracting 'Gene Expression' slot")
        mat <- mat[["Gene Expression"]]
    }
    mat
}

message("Loading raw matrix: ", args$raw_h5)
raw_mat <- load_gex(args$raw_h5)

message("Loading filtered matrix: ", args$filtered_h5)
filt_mat <- load_gex(args$filtered_h5)

message("Loading DropletQC table: ", args$cell_qc)
cell_qc <- read.table(args$cell_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(cell_qc) <- cell_qc$barcode

keep <- cell_qc$barcode[cell_qc$cell_status == "cell"]
if (length(keep) == 0) {
    stop("DropletQC removed all cells; cannot continue with SoupX.")
}

missing_keep <- setdiff(keep, colnames(filt_mat))
if (length(missing_keep) > 0) {
    stop("Some DropletQC barcodes were not found in the filtered matrix: ",
         paste(head(missing_keep, 5), collapse = ", "))
}

filt_mat <- filt_mat[, keep, drop = FALSE]
message("Cells remaining after DropletQC filter: ", ncol(filt_mat))

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

message("Creating SoupX SoupChannel object...")
sc_obj <- SoupChannel(tod = raw_mat, toc = filt_mat)
sc_obj <- setClusters(sc_obj, clusters)

if (!is.null(args$contaminant_genes) && length(args$contaminant_genes) > 0) {
    present_genes <- intersect(args$contaminant_genes, rownames(filt_mat))
    if (length(present_genes) > 0) {
        message("Using user-supplied contaminant genes: ", paste(present_genes, collapse = ", "))
        sc_obj <- setContaminationFraction(sc_obj, useToEst = present_genes)
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

rho_df <- sc_obj$metaData[, "rho", drop = FALSE]
rho_df$barcode <- rownames(rho_df)
colnames(rho_df)[colnames(rho_df) == "rho"] <- "soupx_rho"
cell_qc <- merge(cell_qc, rho_df, by = "barcode", all.x = TRUE)
rownames(cell_qc) <- cell_qc$barcode
cell_qc$soupx_rho[cell_qc$cell_status == "damaged_cell"] <- NA_real_
if (!"nf_umi" %in% colnames(cell_qc)) {
    stop("DropletQC cell_qc.tsv is missing the required nf_umi column.")
}

message("Writing corrected matrix to: ", file.path(args$output_dir, "corrected"))
write10xCounts(
    path = file.path(args$output_dir, "corrected"),
    x = corrected_mat,
    type = "sparse",
    overwrite = TRUE
)

cell_qc_out <- cell_qc[, c("barcode", "soupx_rho", "cell_status", "nf_umi"), drop = FALSE]
cell_qc_path <- file.path(args$output_dir, "cell_qc.tsv")
write.table(cell_qc_out, cell_qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written cell QC to: ", cell_qc_path)

soup_prof <- sc_obj$soupProfile[order(-sc_obj$soupProfile$est), ]
soup_prof$gene <- rownames(soup_prof)
soup_path <- file.path(args$output_dir, "soup_profile.tsv")
write.table(soup_prof[, c("gene", "est", "cnts")], soup_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Written soup profile to: ", soup_path)

message("Done.")
