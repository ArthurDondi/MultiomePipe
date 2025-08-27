import argparse
import timeit
from cellbender.remove_background.downstream import load_anndata_from_input_and_output
import anndata as ad

def initialize_parser():
    parser = argparse.ArgumentParser(description='Merging splitted outputs from cellbender')
    parser.add_argument('--cellbender_input', nargs="+", required=True, help="List of input files of cellbender")
    parser.add_argument('--cellbender_output', nargs="+", required=True, help="List of filtered output files of cellbender")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    cellbender_input = args.cellbender_input
    cellbender_output = args.cellbender_output
    output_file = args.output

    # 1. Merge data
    print("1. Load data")
    start = timeit.default_timer()
    adatas = [load_anndata_from_input_and_output(input_file=f1,
                                                 output_file=f2,
                                                 input_layer_key='raw'
                                                 )
                                                for f1, f2 in zip(cellbender_input, cellbender_output)
                ]

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
