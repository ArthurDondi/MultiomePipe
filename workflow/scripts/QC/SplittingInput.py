import sys
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import argparse
import timeit
import h5py
from scipy.sparse import csr_matrix
from pathlib import Path

def plot_elbow(adata, pngout):
    """
    Generates an elbow plot of total UMIs per cell.
    """
    # Compute QC metrics if not already computed
    if 'total_counts' not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # Sort cells by total UMIs
    sorted_umi = np.sort(adata.obs['total_counts'].values)[::-1]

    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(range(1, len(sorted_umi)+1), sorted_umi)
    plt.xlabel('Cells (ordered by total counts)')
    plt.ylabel('Total UMIs per cell')
    plt.title('Elbow plot of total UMIs')
    plt.xscale('log')
    plt.yscale('log')  
    plt.grid(True, which="both", ls="--", lw=0.5)
    # Save plot
    plt.savefig(pngout, dpi=300)
    plt.close()
    print(f"Elbow plot saved to {pngout}")

def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    
    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """
    
    if '.h5' not in file: file = f'{file}.h5'
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))

def subset(adata, outdir, n_cells):
    """
    Subset a fraction of cells and save in 10X format for CellBender.
    """
    n_cells_tot = adata.n_obs
    print(f"Total cells: {n_cells_tot}")
    # Targetting < {n_cells} cells per chunks
    frac = float(n_cells_tot/n_cells)
    # 1 chunk only:
    if frac<=1:
        write_10X_h5(subset, f"{outdir}/split_0.raw_feature_bc_matrix.h5")
    # Multiple chuncks:
    else:
        n_splits = int(n_cells_tot/n_cells)+1
        # Split indices into n_cells chunks
        indices = np.array_split(np.arange(n_cells), n_splits)

        # Create and save subsets
        for i, idx in enumerate(indices, start=1):
            subset = adata[idx, :].copy()  # subset cells, keep all features
            write_10X_h5(subset, f"{outdir}/split_{i}.raw_feature_bc_matrix.h5")
            print(f"Saved split_{i}.raw_feature_bc_matrix.h5 with {subset.n_obs} cells")    

def main():
    parser = argparse.ArgumentParser(description="Generate elbow plot from an AnnData object")
    parser.add_argument("--h5in", type=str, help="Path to the .h5ad AnnData file")
    parser.add_argument("--outdir", type=str, help="Directory for the splitted h5s")
    parser.add_argument("--plotdir", type=str, help="Directory for plots")
    parser.add_argument("--sample", type=str, help="Sample name")
    parser.add_argument("--n_cells", type=int, help="Number of chunks to divide the sample")

    args = parser.parse_args()

    # Load AnnData
    if args.h5in.split('.')[-1] == 'h5':
        adata = sc.read_10x_h5(args.h5in)
    elif args.h5in.split('.')[-1] == 'h5ad':
        adata = ad.io.read_h5ad(args.h5in)
    else:
        print("Input data must either be raw_feature_bc_matrix.h5 or an AnnData .h5ad")
        sys.exit(0) 

    # Plot elbow
    pngout = args.plotdir + args.sample + ".ElbowPlot.png"
    plot_elbow(adata,pngout)

    # Subset data to run on local machine
    
    subset(adata, args.outdir, args.n_cells)

if __name__ == "__main__":
    start_time = timeit.default_timer()
    main()
    end_time = timeit.default_timer()
    print(f"Execution time: {end_time - start_time:.2f} seconds")
