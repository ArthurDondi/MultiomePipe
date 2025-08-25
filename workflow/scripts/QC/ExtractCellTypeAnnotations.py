import argparse
import timeit
import anndata as ad

def initialize_parser():
    parser = argparse.ArgumentParser(description='Extracting adata.obs')
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    adata = ad.io.read_h5ad(input_file)

    adata.obs.to_csv(output_file, sep=",", index=True)  

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
