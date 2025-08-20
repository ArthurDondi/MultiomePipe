import argparse
import timeit
from cellbender.remove_background.downstream import anndata_from_h5
import anndata as ad

def initialize_parser():
    parser = argparse.ArgumentParser(description='Merging splitted outputs from cellbender')
    parser.add_argument('--input', nargs="+", required=True)
    parser.add_argument('--output', type=str, required=True)
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_files = args.input
    output_file = args.output

    # 1. Merge data
    print("1. Load data")
    start = timeit.default_timer()
    adatas = [anndata_from_h5(f) for f in input_files]
    for adata in adatas:
        adata.var_names_make_unique()
    adata_concat = ad.concat(
        adatas,
        merge="first",
        uns_merge="first",
        join="outer",
        label="split",
        keys=[f"split{i}" for i in range(len(adatas))]
    )
    adata_concat.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
