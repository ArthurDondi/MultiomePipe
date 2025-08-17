#!/usr/bin/env Rscript

suppressMessages(library(SoupX))
suppressMessages(library(scDblFinder))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Matrix))
suppressMessages(library(optparse))

# -----------------------------
# 1. Command-line arguments
# -----------------------------
option_list <- list(
  make_option(c("--raw_dir"), type="character", help="raw 10X matrix folder"),
  make_option(c("--filtered_dir"), type="character", help="Filtered (manually made) 10X matrix folder"),
  make_option(c("--clusters"), type="character", help="Scanpy clusters TSV file"),
  make_option(c("--output"), type="character", help="Output RDS file prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

raw_dir <- opt$raw_dir
filtered_dir <- opt$filtered_dir
cluster_file <- opt$clusters
output_file <- opt$output

# -----------------------------
# 2. Load 10X matrix
# -----------------------------
cat("\nLoading 10X filtered matrix...\n")
# Read matrix, barcodes, features manually
filtered_counts <- readMM(gzfile(file.path(filtered_dir, "matrix.mtx.gz")))
filtered_barcodes <- read.table(gzfile(file.path(filtered_dir, "barcodes.tsv.gz")), stringsAsFactors = FALSE)[,1]
filtered_features <- read.table(gzfile(file.path(filtered_dir, "features.tsv.gz")), stringsAsFactors = FALSE)[,1]
# Set row/column names
rownames(filtered_counts) <- filtered_features
colnames(filtered_counts) <- filtered_barcodes

cat("\nLoading 10X raw matrix (table of droplets)...\n")
# Read matrix, barcodes, features manually
raw_counts <- readMM(gzfile(file.path(raw_dir, "matrix.mtx.gz")))
raw_barcodes <- read.table(gzfile(file.path(raw_dir, "barcodes.tsv.gz")), stringsAsFactors = FALSE)[,1]
raw_features <- read.table(gzfile(file.path(raw_dir, "features.tsv.gz")),       
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE
    )[,1]


# Set row/column names
rownames(raw_counts) <- raw_features
colnames(raw_counts) <- raw_barcodes

# Find indices of filtered_features in raw_features
idx <- match(filtered_features, raw_features)

# Subset raw_counts and raw_features
raw_counts <- raw_counts[idx, , drop = FALSE]
raw_features <- raw_features[idx]

# Sanity check
cat(all(raw_features == filtered_features))  # should be TRUE


#sce <- SingleCellExperiment(assays = list(counts = filtered_counts))
#cat("Cells:", ncol(sce), "Genes:", nrow(sce), "\n")

# -----------------------------
# 3. Load Scanpy clusters
# -----------------------------
cat("\nLoading Scanpy clusters...\n")
scanpy_clusters <- read.table(cluster_file, header=FALSE, stringsAsFactors = FALSE)
colnames(scanpy_clusters) <- "cluster"

# if (nrow(scanpy_clusters) != ncol(sce)) {
#   stop("Number of cells in cluster file does not match SCE object!")
# }

# -----------------------------
# 4. SoupX
# -----------------------------
cat("\nRunning SoupX...\n")
sc <- SoupChannel(tod=raw_counts, toc=filtered_counts)
sc <- setClusters(sc, scanpy_clusters$cluster)
sc <- autoEstCont(sc)

adjusted_counts <- adjustCounts(sc)


assay(sce, "adjusted_counts") <- adjusted_counts

# -----------------------------
# 5. scDblFinder
# -----------------------------
cat("\nRunning scDblFinder...\n")
sce <- scDblFinder(sce)

# -----------------------------
# 6. Save SCE object
# -----------------------------
cat("\nSaving SCE object...\n")
saveRDS(sce, file = output_file)
cat("Done!\n")
