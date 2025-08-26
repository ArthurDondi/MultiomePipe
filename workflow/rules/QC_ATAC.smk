# Snakefile
from pathlib import Path
import os
import pandas as pd

rule all_QC_ATAC:
    input:
        expand("QC/ATAC/{sample}_consensus_peaks.bed", sample=SAMPLES)

rule ExportPseudobulk:
    input:
        fragments = f"{INPUT}/raw/{{sample}}.fragments.tsv.gz",
        cell_data = "QC/RNA/{sample}/Annotation/{sample}_CellTypeAnnotations.csv"
    output:
        bed_paths = "QC/ATAC/{sample}/bed_paths.tsv",
        bw_paths  = "QC/ATAC/{sample}/bw_paths.tsv"
    params:
        script = f"{workflow.basedir}/scripts/ExportPseudobulk.py",
        outdir = lambda wildcards: f"QC/ATAC/{wildcards.sample}",
        variable = config["QC_ATAC"]["ExportPseudobulk"]["variable"],
        sample_id_col = config["QC_ATAC"]["ExportPseudobulk"]["sample_id_col"],
        n_cpu = config["QC_ATAC"]["ExportPseudobulk"]["n_cpu"],
        normalize_bigwig = config["QC_ATAC"]["ExportPseudobulk"]["normalize_bigwig"],
        split_pattern = config["QC_ATAC"]["ExportPseudobulk"]["split_pattern"],
        chromsizes = config["QC_ATAC"]["ExportPseudobulk"]["chromsizes"]
    conda:
        "../envs/pycistopic.yaml"
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
            --outdir {params.outdir} \
            --variable {params.variable} \
            --sample_id_col {params.sample_id_col} \
            --n_cpu {params.n_cpu} \
            {params.normalize_bigwig} \
            --split_pattern {params.split_pattern} \
            {("--chromsizes " + params.chromsizes) if params.chromsizes else ""}
        """


###############################################################################
# Helper: read bed_paths.tsv produced by export_pseudobulk checkpoint
###############################################################################
def _get_celltypes_for_sample(wildcards):
    ck = checkpoints.export_pseudobulk.get(sample=wildcards.sample)
    bed_paths_file = ck.output.bed_paths
    celltypes = []
    with open(bed_paths_file) as fh:
        for line in fh:
            name, path = line.strip().split("\t", 1)
            celltypes.append(name)
    return celltypes

def _get_macs_outputs(wildcards):
    # produce list of expected MACS narrowPeak files (one per celltype)
    celltypes = _get_celltypes_for_sample(wildcards)
    return expand("QC/ATAC/MACS/{{sample}}/{{sample}}_{celltype}_peaks.narrowPeak", celltype=CTYPES[wildcards.sample])

###############################################################################
# 2) MACS2 per sample × celltype (expanded from the consensus rule below)
###############################################################################
rule MACS2CallPeaks:
    input:
        fragments = lambda wildcards: os.path.join(
            "outs/consensus_peak_calling", wildcards.sample,
            "pseudobulk_bed_files", f"{wildcards.celltype}.fragments.tsv.gz"
        )
    output:
        "QC/ATAC/MACS/{sample}/{sample}_{celltype}_peaks.narrowPeak"
    params:
        outdir = lambda wildcards: f"QC/ATAC/MACS/{wildcards.sample}",
        format = config["QC_ATAC"]["MACS2CallPeaks"]["format"],
        gsize = config["QC_ATAC"]["MACS2CallPeaks"]["gsize"],
        qvalue = config["QC_ATAC"]["MACS2CallPeaks"]["qvalue"],
        shift = config["QC_ATAC"]["MACS2CallPeaks"]["shift"],
        extsize = config["QC_ATAC"]["MACS2CallPeaks"]["extsize"],
        keep_dup = config["QC_ATAC"]["MACS2CallPeaks"]["keep_dup"],
        call_summits = config["QC_ATAC"]["MACS2CallPeaks"]["call_summits"],
        nolambda = config["QC_ATAC"]["MACS2CallPeaks"]["nolambda"]
    conda:
        "../envs/macs2.yaml"
    threads: 1
    log:
        "logs/MACS2CallPeaks/{sample}_{celltype}.log"
    benchmark:
        "benchmark/MACS2CallPeaks/{sample}_{celltype}.benchmark.txt"
    shell:
        r"""
        mkdir -p {params.outdir}
        macs2 callpeak \
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
        "QC/ATAC/{sample}_consensus_peaks.bed"
    params:
        script = f"{workflow.basedir}/scripts/GetConsensusPeaks.py"
    conda:
        "../envs/pycistopic.yaml"
    log:
        "logs/GetConsensusPeaks/{sample}.log"
    benchmark:
        "benchmark/GetConsensusPeaks/{sample}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
            --inputs {input.narrowpeaks} \
            --output {output} \
            --sample {wildcards.sample}
        """
