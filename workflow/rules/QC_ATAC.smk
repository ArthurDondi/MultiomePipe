# Snakefile
from pathlib import Path
import os
import pandas as pd

### This rule is only there because I cannot process a whole ATAC-seq dataset locally with 16GB RAM
SUBSAMPLE = True
rule Subsample:
    input:
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
    output:
        subsampled = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz",
        tbi = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz.tbi",
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/Subsample.py",
        fraction = 0.1,
    shell:
        r"""
        python {params.script} \
            --input {input.fragments} \
            --output {output.subsampled} \
            --fraction {params.fraction}
        """

checkpoint ExportPseudobulk:
    input:
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz" if SUBSAMPLE else f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
        cell_data = "QC/RNA/Merged/Annotation/merged.CellTypeAnnotations.csv"
    output:
        bed_paths = "QC/ATAC/{sample}/PseudoBulk/bed_paths.tsv",
        bw_paths  = "QC/ATAC/{sample}/PseudoBulk/bw_paths.tsv"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/ExportPseudoBulk.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/PseudoBulk",
        genome = config['General']['genome'],
        celltype_key = config['General']['celltype_key'],
        sample_key = config['General']["sample_key"],
        n_cpu = config["QC_ATAC"]["ExportPseudobulk"]["n_cpu"],
        normalize_bigwig = config["QC_ATAC"]["ExportPseudobulk"]["normalize_bigwig"],
        split_pattern = config["QC_ATAC"]["ExportPseudobulk"]["split_pattern"],
        chromsizes = lambda wildcards: f"--chromsizes {config['QC_ATAC']['ExportPseudobulk']['chromsizes']}" 
                                       if config['QC_ATAC']['ExportPseudobulk']['chromsizes'] else "",
    threads: config["QC_ATAC"]["ExportPseudobulk"]["n_cpu"]
    log:
        "logs/ExportPseudobulk/{sample}.log"
    benchmark:
        "benchmark/ExportPseudobulk/{sample}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
            --fragments {input.fragments} \
            --cell_data {input.cell_data} \
            --sample {wildcards.sample} \
            --genome {params.genome} \
            --outdir {params.outdir} \
            --celltype_key {params.celltype_key} \
            --sample_key {params.sample_key} \
            --n_cpu {params.n_cpu} \
            --split_pattern {params.split_pattern} {params.normalize_bigwig} {params.chromsizes}
        """

###############################################################################
# Helper: read bed_paths.tsv produced by export_pseudobulk checkpoint
###############################################################################
def _get_celltypes_for_sample(wildcards):
    ck = checkpoints.ExportPseudobulk.get(sample=wildcards.sample)
    bed_paths_file = ck.output.bed_paths
    celltypes = []
    with open(bed_paths_file) as fh:
        for line in fh:
            name, path = line.strip().split("\t", 1)
            celltypes.append(name.replace(' ','_'))
    return celltypes

def _get_macs_outputs(wildcards):
    # produce list of expected MACS narrowPeak files (one per celltype)
    CTYPES = _get_celltypes_for_sample(wildcards)
    return expand("QC/ATAC/{sample}/MACS/{sample}_{celltype}_peaks.narrowPeak",
                    sample = wildcards.sample,
                    celltype=CTYPES)

###############################################################################
# 2) MACS2 per sample × celltype (expanded from the consensus rule below)
###############################################################################
rule MACS2CallPeaks:
    input:
        fragments = "QC/ATAC/{sample}/PseudoBulk/pseudobulk_bed_files/{celltype}.fragments.tsv.gz"
    output:
        narrowpeak = "QC/ATAC/{sample}/MACS/{sample}_{celltype}_peaks.narrowPeak"
    params:
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/MACS",
        format = config["QC_ATAC"]["MACS2CallPeaks"]["format"],
        gsize = config["QC_ATAC"]["MACS2CallPeaks"]["gsize"],
        qvalue = config["QC_ATAC"]["MACS2CallPeaks"]["qvalue"],
        shift = config["QC_ATAC"]["MACS2CallPeaks"]["shift"],
        extsize = config["QC_ATAC"]["MACS2CallPeaks"]["extsize"],
        keep_dup = config["QC_ATAC"]["MACS2CallPeaks"]["keep_dup"],
        call_summits = config["QC_ATAC"]["MACS2CallPeaks"]["call_summits"],
        nolambda = config["QC_ATAC"]["MACS2CallPeaks"]["nolambda"]
    threads: 1
    log:
        "logs/MACS2CallPeaks/{sample}_{celltype}.log"
    benchmark:
        "benchmark/MACS2CallPeaks/{sample}_{celltype}.benchmark.txt"
    shell:
        r"""
        mkdir -p {params.outdir}
        macs3 callpeak \
            --treatment {input.fragments} \
            --name {wildcards.sample}_{wildcards.celltype} \
            --outdir {params.outdir} \
            --format {params.format} \
            --gsize {params.gsize} \
            --qvalue {params.qvalue} \
            --nomodel \
            --shift {params.shift} \
            --extsize {params.extsize} \
            --keep-dup {params.keep_dup} \
            {params.call_summits} \
            {params.nolambda} 
        """

###############################################################################
# 3) Consensus peaks: wait for all MACS outputs, then run pycisTopic merging
###############################################################################
rule GetConsensusPeaks:
    input:
        narrowpeaks = lambda wildcards: _get_macs_outputs(wildcards)
    output:
        consensus = "QC/ATAC/{sample}/ConsensusPeaks/{sample}_consensus_peaks.bed"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/GetConsensusPeaks.py",
        genome = config['General']['genome'],
        path_to_blacklist = config["QC_ATAC"]["GetConsensusPeaks"]["path_to_blacklist"]
    log:
        "logs/GetConsensusPeaks/{sample}.log"
    benchmark:
        "benchmark/GetConsensusPeaks/{sample}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
            --inputs {input.narrowpeaks} \
            --output {output} \
            --sample {wildcards.sample} \
            --genome {params.genome} \
            --path_to_blacklist {params.path_to_blacklist} 
        """

rule DownloadTSS:
    output:
        tss_annotation = "QC/ATAC/TSS/tss_annotation.bed"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/DownloadTSS.py",
        tss_dataset = config["QC_ATAC"]["DownloadTSS"]['tss_dataset'],
        tss_host = config["QC_ATAC"]["DownloadTSS"]['tss_host'],
    log:
        "logs/DownloadTSS/AllSamples.log"
    benchmark:
        "benchmark/DownloadTSS/AllSamples.benchmark.txt"
    shell:
        r"""
        python {params.script} \
            --tss_dataset {params.tss_dataset} \
            --tss_host {params.tss_host} \
            --output {output.tss_annotation}
        """

rule ComputeQC:
    input:
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz" if SUBSAMPLE else f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
        consensus = "QC/ATAC/{sample}/ConsensusPeaks/{sample}_consensus_peaks.bed",
        tss_annotation = "QC/ATAC/TSS/tss_annotation.bed",
    output:
        parquet = "QC/ATAC/{sample}/QC/{sample}.fragments_insert_size_dist.parquet",
    params:
        out_prefix = lambda wildcards: f"QC/ATAC/{wildcards.sample}/QC/{wildcards.sample}",
        min_fragments_per_cb = config["QC_ATAC"]['ComputeQC']['min_fragments_per_cb'],
        dont_collapse_duplicates = config["QC_ATAC"]['ComputeQC']['dont-collapse_duplicates'],
        tss_flank_window = config["QC_ATAC"]['ComputeQC']['tss_flank_window'],
        tss_smoothing_rolling_window = config["QC_ATAC"]['ComputeQC']['tss_smoothing_rolling_window'],
        tss_minimum_signal_window = config["QC_ATAC"]['ComputeQC']['tss_minimum_signal_window'],
        tss_window = config["QC_ATAC"]['ComputeQC']['tss_window'],
        tss_min_norm = config["QC_ATAC"]['ComputeQC']['tss_min_norm'],
    threads: config["QC_ATAC"]['ComputeQC']['threads']
    log:
        "logs/ComputeQC/{sample}.log"
    benchmark:
        "benchmark/ComputeQC/{sample}.benchmark.txt"
    shell:
        r"""
        pycistopic qc \
            --fragments {input.fragments} \
            --regions {input.consensus} \
            --tss {input.tss_annotation} \
            --output {params.out_prefix} \
            --threads {threads} \
            --min_fragments_per_cb {params.min_fragments_per_cb} \
            --tss_flank_window {params.tss_flank_window} \
            --tss_smoothing_rolling_window {params.tss_smoothing_rolling_window} \
            --tss_minimum_signal_window {params.tss_minimum_signal_window} \
            --tss_window {params.tss_window} \
            --tss_min_norm {params.tss_min_norm} {params.dont_collapse_duplicates}
        """

rule ApplyQC:
    input:
        parquet = "QC/ATAC/{sample}/QC/{sample}.fragments_insert_size_dist.parquet",
    output:
        barcodes = "QC/ATAC/{sample}/QC/{sample}.barcodes_passing_qc.tsv",
        png ="QC/ATAC/{sample}/QC/Plots/sample_statistics.png"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/ApplyQC.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/QC",
        unique_fragments_threshold = config["QC_ATAC"]["ApplyQC"]['unique_fragments_threshold'],
        tss_enrichment_threshold = config["QC_ATAC"]["ApplyQC"]['tss_enrichment_threshold'],
        frip_threshold = config["QC_ATAC"]["ApplyQC"]['frip_threshold'],
    shell:
        r"""
        python {params.script} \
            --outdir {params.outdir} \
            --sample {wildcards.sample} \
            --barcodes {output.barcodes} \
            --unique_fragments_threshold {params.unique_fragments_threshold} \
            --tss_enrichment_threshold {params.tss_enrichment_threshold} \
            --frip_threshold {params.frip_threshold}
        """

