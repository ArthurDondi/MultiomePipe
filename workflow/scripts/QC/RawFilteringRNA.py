import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import numpy as np
import json

# -----------------------------
# Functions
# -----------------------------
def load_data(input_file,sample,donor):
    adata = ad.io.read_h5ad(input_file)
    # Adding sample and donor metadata
    adata.obs['sample'] = sample
    adata.obs['donor'] = donor
    if "cellbender" in adata.layers: 
        adata.X = adata.layers['cellbender']
    elif "counts" in adata.layers: 
        adata.X = adata.layers['counts']
    else:
        adata.layers["counts"] = adata.X.copy()
    adata.var_names_make_unique()
    return adata

def run_calculate_qc_metrics(adata, sample):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=["mt","ribo","hb"],
                               inplace=True, 
                               log1p=True)
    adata.var_names_make_unique()

    # Plotting
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_{sample}.png"
    )
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_{sample}.png"
    )
    return adata

def run_raw_filtering(adata, min_genes, min_cells, n_mads, sample):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    #Filtering cells based on mitochondrial reads %, by n_mads Median Absolute Deviation (MAD)
    x = adata.obs['pct_counts_mt']
    median = np.median(x)
    mad = np.median(np.abs(x - median))
    mask = np.abs(x - median) <= n_mads * mad
    adata= adata[mask].copy()

    print(f"Mito % Median: {median:.2f}, MAD: {mad:.2f}")
    print(f"Cells kept: {mask.sum()} / {len(mask)}, ({100 * mask.sum() / len(mask):.2f}%)")

    # Plotting# Plotting
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_{min_genes}minGenes_{min_cells}minCells_{sample}.png"
    )
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_{min_genes}minGenes_{min_cells}minCells_{sample}.png"
    )
    return adata

def identify_doublets(adata, doublet_threshold, sample):
    total_cells = adata.n_obs
    sc.pp.scrublet(adata, doublet_threshold)
    
    n_predicted_doublets = adata.obs['predicted_doublet'].sum()
    
    print(f"Predicted doublets: {n_predicted_doublets} / {total_cells} cells ({100 * n_predicted_doublets / total_cells:.1f}%)")

    sc.pl.scrublet_score_distribution(adata,
                                      show=False,
                                      save=f"_{sample}.png")

    # Filtering doublets
    adata = adata[~adata.obs["predicted_doublet"]]

    return adata


def run_normalization_and_clustering(adata, sample, n_top_genes, is_filtered):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # Define cell cycle gene lists (human genes)
    s_genes = [
        'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6',
        'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1',
        'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2',
        'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1',
        'CHAF1B', 'BRIP1', 'E2F8'
    ]

    g2m_genes = [
        'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2',
        'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2',
        'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B',
        'HJURP', 'CDCA3', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5',
        'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5',
        'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
    ]

    # Score cell cycle phases
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    # Plot highly variable genes:
    sc.pl.highly_variable_genes(adata,
                                show=False,
                                save=f"_HVG_{sample}.png")

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    # Plot PCs variance ratio :
    sc.pl.pca_variance_ratio(adata, 
                             n_pcs=50, 
                             log=True,
                             show=False,
                             save=f"_{sample}.png")
    sc.tl.umap(adata)
    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
    sc.pl.umap(adata,
            color=[f"leiden_res_{res:4.2f}" for res in [0.2, 0.5, 1]],
            legend_loc="on data",
            show=False,
            save=f"_leiden_res_{sample}.png")
    
    
    QCs = ["log1p_total_counts", "n_genes_by_counts", "pct_counts_mt", "doublet_score",'S_score','G2M_score']
    if not is_filtered:
        QCs.append("background_fraction")
    sc.pl.umap(adata,
                color=QCs,
                ncols=3,
                show=False,
                save=f"_QC_{sample}.png")

    return adata


def initialize_parser():
    parser = argparse.ArgumentParser(description='QC + clustering for 10X Multiome RNA')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument('--donor', type=str, required=True)
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--min_genes', type=int, default=100)
    parser.add_argument('--min_cells', type=int, default=3)
    parser.add_argument('--n_mads', type=int, default=3)
    parser.add_argument('--doublet_threshold', type=float, default=.3)
    parser.add_argument('--n_top_genes', type=int, default=2000)
    parser.add_argument('--is_filtered', action='store_true')
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    sample = args.sample
    donor = args.donor
    plotdir = args.plotdir
    min_genes = args.min_genes
    min_cells = args.min_cells
    n_mads = args.n_mads
    doublet_threshold = args.doublet_threshold
    n_top_genes = args.n_top_genes
    is_filtered = args.is_filtered

    sc.settings.figdir = plotdir

    # 1. Load data
    print("1. Load data")
    start = timeit.default_timer()
    adata = load_data(input_file,sample,donor)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

    # 2. QC metrics
    print("2. QC metrics")
    start = timeit.default_timer()
    adata = run_calculate_qc_metrics(adata, sample)
    stop = timeit.default_timer()
    print(f"QC metrics computed in {round(stop-start,2)}s")

    # 3. RaW Filtering
    print("3. Filtering")
    start = timeit.default_timer()
    adata = run_raw_filtering(adata, min_genes, min_cells, n_mads, sample)
    stop = timeit.default_timer()
    print(f"Filtering done in {round(stop-start,2)}s")

    # 4. Detecting Doublets
    print("4. Detecting Doublets")
    start = timeit.default_timer()
    adata = identify_doublets(adata,doublet_threshold, sample)
    stop = timeit.default_timer()
    print(f"Detecting doublets done in {round(stop-start,2)}s")
    
    # 5. Normalization + clustering
    print("5. Normalization + clustering")
    start = timeit.default_timer()
    adata = run_normalization_and_clustering(adata, sample, n_top_genes, is_filtered)
    stop = timeit.default_timer()
    print(f"Normalization + clustering done in {round(stop-start,2)}s")

    # 7. Write data
    print("7. Write data")
    start = timeit.default_timer()
    adata.write(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
