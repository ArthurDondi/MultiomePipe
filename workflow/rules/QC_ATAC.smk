# Snakefile
from pathlib import Path
import os
import pandas as pd

### This rule is only there because I cannot process a whole ATAC-seq dataset locally with 16GB RAM
SUBSAMPLE = True
rule Subsample:
    input:
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
        cell_data = "QC/RNA/Merged/Annotation/merged.CellTypeAnnotations.csv"
    output:
        subsampled = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz",
        tbi = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz.tbi",
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/Subsample.py",
        fraction = 0.5,
    shell:
        r"""
        python {params.script} \
            --input {input.fragments} \
            --output {output.subsampled} \
            --cell_data {input.cell_data} \
            --sample {wildcards.sample} \
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
        exec > {log} 2>&1
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

rule MACS3CallPeaks:
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
        "logs/MACS3CallPeaks/{sample}_{celltype}.log"
    benchmark:
        "benchmark/MACS3CallPeaks/{sample}_{celltype}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
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

rule GetConsensusPeaks:
    input:
        narrowpeaks = lambda wildcards: _get_macs_outputs(wildcards)
    output:
        consensus = "QC/ATAC/{sample}/ConsensusPeaks/{sample}_consensus_peaks.bed"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/GetConsensusPeaks.py",
        genome = config['General']['genome'],
        path_to_blacklist = f"{workflow.basedir}/../data/hg38-blacklist.v2.bed",
    log:
        "logs/GetConsensusPeaks/{sample}.log"
    benchmark:
        "benchmark/GetConsensusPeaks/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
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
        exec > {log} 2>&1
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
        parquet = "QC/ATAC/{sample}/QC/{sample}.fragments_stats_per_cb.parquet",
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
        exec > {log} 2>&1
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
        parquet = "QC/ATAC/{sample}/QC/{sample}.fragments_stats_per_cb.parquet",
    output:
        barcodes = "QC/ATAC/{sample}/QC/{sample}.barcodes_passing_qc.tsv",
        png ="QC/ATAC/{sample}/QC/Plots/sample_statistics.png"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/ApplyQC.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/QC",
        unique_fragments_threshold = config["QC_ATAC"]["ApplyQC"]['unique_fragments_threshold'],
        tss_enrichment_threshold = config["QC_ATAC"]["ApplyQC"]['tss_enrichment_threshold'],
        frip_threshold = config["QC_ATAC"]["ApplyQC"]['frip_threshold'],
    log:
        "logs/ApplyQC/{sample}.log"
    benchmark:
        "benchmark/ApplyQC/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python {params.script} \
            --outdir {params.outdir} \
            --sample {wildcards.sample} \
            --barcodes {output.barcodes} \
            --unique_fragments_threshold {params.unique_fragments_threshold} \
            --tss_enrichment_threshold {params.tss_enrichment_threshold} \
            --frip_threshold {params.frip_threshold}
        """

rule CreateATACCountMatrix:
    input:
        barcodes = "QC/ATAC/{sample}/QC/{sample}.barcodes_passing_qc.tsv",
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.subsample.tsv.gz" if SUBSAMPLE else f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
        consensus = "QC/ATAC/{sample}/ConsensusPeaks/{sample}_consensus_peaks.bed",
        parquet = "QC/ATAC/{sample}/QC/{sample}.fragments_stats_per_cb.parquet",
    output:
        pkl = "QC/ATAC/{sample}/QC/{sample}.matrix.pkl"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/CreateATACCountMatrix.py",
        blacklist = f"{workflow.basedir}/../data/hg38-blacklist.v2.bed",
        split_pattern = config["QC_ATAC"]["ExportPseudobulk"]["split_pattern"],
        min_fragments_per_region = config["QC_ATAC"]['CreateATACCountMatrix']['min_fragments_per_region'],
        min_cells_per_region = config["QC_ATAC"]['CreateATACCountMatrix']['min_cells_per_region'],
    threads: config["QC_ATAC"]['CreateATACCountMatrix']['threads']
    log:
        "logs/CreateATACCountMatrix/{sample}.log"
    benchmark:
        "benchmark/CreateATACCountMatrix/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python {params.script} \
            --sample {wildcards.sample} \
            --barcodes {input.barcodes} \
            --fragments {input.fragments} \
            --regions {input.consensus} \
            --out_cistopic_obj {output.pkl} \
            --blacklist {params.blacklist} \
            --parquet {input.parquet} \
            --split_pattern {params.split_pattern} \
            --min_fragments_per_region {params.min_fragments_per_region} \
            --min_cells_per_region {params.min_cells_per_region} \
            --n_cpu {threads} 
        """

rule RunLDAModels:
    input:
        pkl = "QC/ATAC/{sample}/QC/{sample}.matrix.pkl"
    output:
        pkl = "QC/ATAC/{sample}/LDAModels/{sample}.models.pkl",
        png = "QC/ATAC/{sample}/LDAModels/Plots/model_evaluation.png"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/RunLDAModels.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/LDAModels",
        n_topics = config["QC_ATAC"]['RunLDAModels']['n_topics'],
        n_iter = config["QC_ATAC"]['RunLDAModels']['n_iter'],
        random_state = config["QC_ATAC"]['RunLDAModels']['random_state'],
        alpha = config["QC_ATAC"]['RunLDAModels']['alpha'],
        alpha_by_topic = config["QC_ATAC"]['RunLDAModels']['alpha_by_topic'],
        eta = config["QC_ATAC"]['RunLDAModels']['eta'],
        eta_by_topic = config["QC_ATAC"]['RunLDAModels']['eta_by_topic'],
    threads: config["QC_ATAC"]['RunLDAModels']['n_cpus']
    log:
        "logs/RunLDAModels/{sample}.log"
    benchmark:
        "benchmark/RunLDAModels/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python {params.script} \
            --outdir {params.outdir} \
            --in_cistopic_obj {input.pkl} \
            --out_cistopic_obj {output.pkl} \
            --n_cpu {threads} \
            --n_topics {params.n_topics} \
            --n_iter {params.n_iter} \
            --random_state {params.random_state} \
            --alpha {params.alpha} {params.alpha_by_topic} \
            --eta {params.eta} {params.eta_by_topic} 
        """

rule Clustering:
    input:
        pkl = "QC/ATAC/{sample}/LDAModels/{sample}.models.pkl",
        cell_data = "QC/RNA/Merged/Annotation/merged.CellTypeAnnotations.csv"
    output:
        pkl = "QC/ATAC/{sample}/Clustering/{sample}.clustering.pkl",
        png = "QC/ATAC/{sample}/Clustering/Plots/heatmap.png"
    params:
        script = f"{workflow.basedir}/scripts/QC_ATAC/Clustering.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}/Clustering",
        sample_key = config['General']["sample_key"],
    log:
        "logs/Clustering/{sample}.log"
    benchmark:
        "benchmark/Clustering/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python {params.script} \
            --outdir {params.outdir} \
            --in_cistopic_obj {input.pkl} \
            --out_cistopic_obj {output.pkl} \
            --cell_data {input.cell_data} \
            --sample {wildcards.sample} \
            --sample_key {params.sample_key}
        """
    