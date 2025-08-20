#!/usr/bin/env python3
"""
Export pseudobulk fragments and bigwigs per celltype using pycisTopic.export_pseudobulk

Writes:
 - outs/.../bed_paths.tsv  (each line: <group>\t<path_to_fragments.tsv.gz>)
 - outs/.../bw_paths.tsv   (each line: <group>\t<path_to_bigwig.bw>)
"""
import argparse
import os
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

def load_chromsizes(path=None):
    if path:
        chromsizes = pd.read_table(path, header=None, names=["Chromosome", "End"])
    else:
        url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
        chromsizes = pd.read_table(url, header=None, names=["Chromosome", "End"])
    chromsizes.insert(1, "Start", 0)
    return chromsizes

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fragments", required=True, help="path to fragments.tsv.gz for this sample")
    p.add_argument("--cell_data", required=True, help="cell_data.tsv (index is barcode)")
    p.add_argument("--sample", required=True)
    p.add_argument("--outdir", required=True, help="output dir (per sample)")
    p.add_argument("--variable", default="VSN_cell_type", help="column in cell_data with cell type labels")
    p.add_argument("--sample_id_col", default="VSN_sample_id", help="sample id column in cell_data")
    p.add_argument("--n_cpu", type=int, default=4)
    p.add_argument("--normalize_bigwig", action="store_true")
    p.add_argument("--split_pattern", default="-")
    p.add_argument("--chromsizes", default=None, help="optional chromsizes file (3rd party). If not provided, downloads hg38 chromsizes.")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    bed_path = os.path.join(args.outdir, "pseudobulk_bed_files")
    bw_path  = os.path.join(args.outdir, "pseudobulk_bw_files")
    os.makedirs(bed_path, exist_ok=True)
    os.makedirs(bw_path, exist_ok=True)

    # read cell metadata
    cell_data = pd.read_table(args.cell_data, index_col=0)

    # fragments mapping - pycisTopic expects a dict of sample->path
    fragments_dict = { args.sample: args.fragments }

    chromsizes = load_chromsizes(args.chromsizes)

    # call export_pseudobulk (this mirrors the code in the notebook)
    bw_paths, bed_paths = export_pseudobulk(
        input_data = cell_data,
        variable = args.variable,
        sample_id_col = args.sample_id_col,
        chromsizes = chromsizes,
        bed_path = bed_path,
        bigwig_path = bw_path,
        path_to_fragments = fragments_dict,
        n_cpu = args.n_cpu,
        normalize_bigwig = args.normalize_bigwig,
        temp_dir = "/tmp",
        split_pattern = args.split_pattern
    )

    # write mapping files for downstream rules
    bed_paths_file = os.path.join(args.outdir, "bed_paths.tsv")
    with open(bed_paths_file, "wt") as fh:
        for k, v in bed_paths.items():
            fh.write(f"{k}\t{v}\n")

    bw_paths_file = os.path.join(args.outdir, "bw_paths.tsv")
    with open(bw_paths_file, "wt") as fh:
        for k, v in bw_paths.items():
            fh.write(f"{k}\t{v}\n")

    print(f"Wrote bed_paths -> {bed_paths_file}")
    print(f"Wrote bw_paths  -> {bw_paths_file}")


if __name__ == "__main__":
    main()
