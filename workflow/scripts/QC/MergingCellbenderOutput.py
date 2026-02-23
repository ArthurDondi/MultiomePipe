import argparse
import timeit
from cellbender.remove_background.downstream import load_anndata_from_input_and_output

def initialize_parser():
    parser = argparse.ArgumentParser(description='Convert cellbender output to h5ad')
    parser.add_argument('--cellbender_input', type=str, required=True, help="Input file of cellbender (raw h5)")
    parser.add_argument('--cellbender_output', type=str, required=True, help="Filtered output file of cellbender (h5)")
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

    # 1. Load data
    print("1. Load data")
    start = timeit.default_timer()
    adata = load_anndata_from_input_and_output(input_file=cellbender_input,
                                               output_file=cellbender_output,
                                               input_layer_key='raw'
                                               )
    adata.var_names_make_unique()
    adata.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
