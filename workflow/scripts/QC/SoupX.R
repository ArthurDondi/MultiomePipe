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
    library(Seurat)
    library(ggplot2)
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

# Comprehensive references (e.g. GENCODE) have non-unique gene symbols;
# Read10X_h5(use.names = TRUE) uniquifies them for processing. Capture the
# Ensembl gene IDs (always unique, read in the same feature order) so the
# corrected matrix can be written with them — downstream alignment keys on IDs.
load_gene_ids <- function(h5_path) {
    ids <- Read10X_h5(h5_path, use.names = FALSE)
    if (is.list(ids)) ids <- ids[["Gene Expression"]]
    rownames(ids)
}
sym_to_id <- setNames(load_gene_ids(args$filtered_h5), rownames(filt_mat))

message("Loading DropletQC table: ", args$cell_qc)
cell_qc <- read.table(args$cell_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(cell_qc) <- cell_qc$barcode
if (!"nf_umi" %in% colnames(cell_qc)) {
    stop("DropletQC cell_qc.tsv is missing the required nf_umi column.")
}

keep <- cell_qc$barcode[cell_qc$cell_status == "cell"]
if (length(keep) == 0) {
    stop("DropletQC removed all cells; cannot continue with SoupX.")
}

missing_keep <- setdiff(keep, colnames(filt_mat))
if (length(missing_keep) > 0) {
    warning(sprintf(
        "%d DropletQC 'cell' barcodes not present in the Cell Ranger filtered matrix and will be skipped: %s%s",
        length(missing_keep),
        paste(head(missing_keep, 5), collapse = ", "),
        if (length(missing_keep) > 5) " ..." else ""
    ))
    keep <- intersect(keep, colnames(filt_mat))
    if (length(keep) == 0) {
        stop("No DropletQC 'cell' barcodes overlap the filtered matrix; cannot continue.")
    }
}

filt_mat <- filt_mat[, keep, drop = FALSE]
message("Cells remaining after DropletQC filter: ", ncol(filt_mat))

message("Running quick Seurat clustering for SoupX...")

# Do NOT apply a min.features cell filter here. The SoupChannel below is built
# from the full DropletQC-passing filt_mat, and SoupX::setClusters requires a
# cluster label for every cell in the channel. Dropping low-gene cells from the
# Seurat object would leave them unlabelled and abort setClusters with the
# opaque "Invalid cluster specification. See help." error. min.cells filters
# genes (not cells), so it is safe and keeps the quick clustering clean.
so <- CreateSeuratObject(counts = filt_mat, min.cells = 3)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = args$resolution, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

clusters <- setNames(as.character(so$seurat_clusters), colnames(so))
message(sprintf("Found %d clusters (resolution = %.1f)", length(unique(clusters)), args$resolution))

# Align cluster labels to the SoupChannel's cells (filt_mat order). Without a
# min.features filter above these sets already match; guard anyway so any future
# regression fails here with a clear message rather than SoupX's opaque one.
unlabelled <- setdiff(colnames(filt_mat), names(clusters))
if (length(unlabelled) > 0) {
    stop(sprintf(
        "%d filtered cell(s) have no Seurat cluster label (e.g. %s); cannot run SoupX::setClusters.",
        length(unlabelled), paste(head(unlabelled, 3), collapse = ", ")
    ))
}
clusters <- clusters[colnames(filt_mat)]

message("Creating SoupX SoupChannel object...")
sc_obj <- SoupChannel(tod = raw_mat, toc = filt_mat)
sc_obj <- setClusters(sc_obj, clusters)
umap_coords <- Embeddings(so, "umap")
sc_obj <- setDR(sc_obj, umap_coords[colnames(filt_mat), , drop = FALSE])

plot_contaminant_genes <- NULL
if (!is.null(args$contaminant_genes) && length(args$contaminant_genes) > 0) {
    plot_contaminant_genes <- intersect(args$contaminant_genes, rownames(filt_mat))
    if (length(plot_contaminant_genes) > 0) {
        message("Contaminant genes for plotChangeMap: ", paste(plot_contaminant_genes, collapse = ", "))
    } else {
        warning("None of the provided contaminant genes found in matrix")
    }
}

message("Estimating contamination with autoEstCont...")
sc_obj <- tryCatch({
    sc_obj <- autoEstCont(sc_obj, verbose = FALSE)
    rho_auto <- median(sc_obj$metaData$rho, na.rm = TRUE)
    message(sprintf("rho (autoEstCont): %.4f (%.2f%%)", rho_auto, rho_auto * 100))
    sc_obj
}, error = function(e) {
    message("Contamination estimation failed (", conditionMessage(e), "), falling back to fixed rho = 0.05")
    setContaminationFraction(sc_obj, 0.05)
})

message("Adjusting counts...")
corrected_mat <- adjustCounts(sc_obj, roundToInt = TRUE, verbose = FALSE)
message(sprintf("SoupX corrected %d genes x %d cells", nrow(corrected_mat), ncol(corrected_mat)))

rho_col <- intersect(c("rho", "contFrac"), colnames(sc_obj$metaData))[1]
if (is.na(rho_col)) {
    warning("Could not find rho/contFrac column in SoupX metaData; soupx_rho will be NA")
    cell_qc$soupx_rho <- NA_real_
} else {
    rho_df <- sc_obj$metaData[, rho_col, drop = FALSE]
    rho_df$barcode <- rownames(rho_df)
    colnames(rho_df)[1] <- "soupx_rho"
    cell_qc <- merge(cell_qc, rho_df, by = "barcode", all.x = TRUE)
    rownames(cell_qc) <- cell_qc$barcode
    cell_qc$soupx_rho[cell_qc$cell_status == "damaged_cell"] <- NA_real_
}

message("Generating plotChangeMap QC plots...")
plot_dir <- file.path(args$output_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

plot_genes <- unique(c("MALAT1", plot_contaminant_genes))
plot_genes <- intersect(plot_genes, rownames(corrected_mat))

for (gene in plot_genes) {
    tryCatch({
        p <- plotChangeMap(sc_obj, corrected_mat, gene) +
            ggtitle(sprintf("%s – SoupX correction", gene))
        ggsave(file.path(plot_dir, sprintf("changeMap_%s.pdf", gene)),
               plot = p, width = 6, height = 5)
    }, error = function(e) {
        warning(sprintf("plotChangeMap failed for %s: %s", gene, conditionMessage(e)))
    })
}
message(sprintf("Written plotChangeMap PDF(s) to: %s", plot_dir))

message("Writing corrected matrix to: ", file.path(args$output_dir, "corrected"))
# Write real Ensembl IDs (unique) as the feature IDs alongside the symbols, so
# the h5ad step can align the raw-count layer on IDs instead of non-unique names.
corrected_ids <- unname(sym_to_id[rownames(corrected_mat)])
if (anyNA(corrected_ids)) {
    stop(sprintf("Could not map %d corrected gene symbol(s) back to Ensembl IDs",
                 sum(is.na(corrected_ids))))
}
write10xCounts(
    path = file.path(args$output_dir, "corrected"),
    x = corrected_mat,
    gene.id = corrected_ids,
    gene.symbol = rownames(corrected_mat),
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
write.table(soup_prof[, c("gene", "est", "counts")], soup_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Written soup profile to: ", soup_path)

message("Done.")