#!/usr/bin/env python3
"""
Collect narrowPeak files (from MACS2), convert each to a PyRanges and call
pycisTopic.iterative_peak_calling.get_consensus_peaks(...).to_bed(output)

Usage:
python GetConsensusPeaks.py --inputs sample_A_celltype1_peaks.narrowPeak sample_A_celltype2_peaks.narrowPeak \
    --output outs/consensus_peak_calling/sampleA_consensus_peaks.bed --sample sampleA
"""
import argparse
import os
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks

def load_chromsizes(path=None):
    if path:
        cs = pd.read_table(path, header=None, names=["Chromosome", "End"])
    else:
        url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
        cs = pd.read_table(url, header=None, names=["Chromosome", "End"])
    cs.insert(1, "Start", 0)
    return cs

def narrowpeak_to_pyranges(path):
    df = pd.read_table(path, header=None)
    # ensure we have column names; narrowPeak spec up to 10 columns
    ncols = df.shape[1]
    base_cols = ["Chromosome", "Start", "End"]
    other_cols = [f"col{i}" for i in range(4, ncols+1)]
    df.columns = base_cols + other_cols[:max(0, ncols-3)]
    # convert to PyRanges
    pr_obj = pr.PyRanges(df)
    return pr_obj, df

def derive_name_from_filename(fname, sample_name=None):
    bn = os.path.basename(fname)
    bn = bn.replace(".narrowPeak", "").replace(".bed", "")
    # if we used naming sample_celltype in Snakemake, remove sample prefix
    if sample_name and bn.startswith(f"{sample_name}_"):
        return bn[len(sample_name)+1:]
    # strip trailing "_peaks" if present
    if bn.endswith("_peaks"):
        bn = bn[:-6]
    return bn

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--sample", required=False, help="optional sample prefix to strip from filenames when naming celltypes")
    p.add_argument("--peak_half_width", type=int, default=250)
    p.add_argument("--chromsizes", default=None)
    p.add_argument("--path_to_blacklist", default=None)
    args = p.parse_args()

    narrow_peaks_dict = {}
    for f in args.inputs:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input narrowPeak not found: {f}")
        pr_obj, df = narrowpeak_to_pyranges(f)
        name = derive_name_from_filename(f, sample_name=args.sample)
        narrow_peaks_dict[name] = pr_obj

    chromsizes = load_chromsizes(args.chromsizes)

    consensus = get_consensus_peaks(
        narrow_peaks_dict,
        peak_half_width = args.peak_half_width,
        chromsizes = chromsizes,
        path_to_blacklist = args.path_to_blacklist
    )

    # write consensus peaks to BED (mirror notebook call)
    outdir = os.path.dirname(args.output)
    os.makedirs(outdir, exist_ok=True)
    consensus.to_bed(path=args.output, keep=True, compression='infer', chain=False)
    print(f"Wrote consensus peaks to {args.output}")

if __name__ == "__main__":
    main()
