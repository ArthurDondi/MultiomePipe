### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets

import os
import json
import warnings

# Background correction method: "cellbender" (default) or "soupx" (DropletQC + SoupX)
BG_CORRECTION = config.get("QC_RNA", {}).get("background_correction", "cellbender")
if BG_CORRECTION == "soupx_dropletqc":
    warnings.warn(
        "QC_RNA.background_correction='soupx_dropletqc' is deprecated; use 'soupx' instead.",
        DeprecationWarning,
        stacklevel=1
    )
    BG_CORRECTION = "soupx"
elif BG_CORRECTION not in {"cellbender", "soupx"}:
    raise ValueError(
        f"Unsupported QC_RNA.background_correction '{BG_CORRECTION}'. "
        "Choose 'cellbender' or 'soupx'."
    )

LEGACY_SOUPX_CONFIG = config.get("QC_RNA", {}).get("SoupXDropletQC")
if LEGACY_SOUPX_CONFIG:
    warnings.warn(
        "QC_RNA.SoupXDropletQC is deprecated; use QC_RNA.DropletQC and QC_RNA.SoupX instead.",
        DeprecationWarning,
        stacklevel=1
    )

SOUPX_CONFIG = (
    config.get("QC_RNA", {}).get("SoupX")
    or LEGACY_SOUPX_CONFIG
    or {}
)
DROPLETQC_CONFIG = (
    config.get("QC_RNA", {}).get("DropletQC")
    or LEGACY_SOUPX_CONFIG
    or {}
)

if BG_CORRECTION == "soupx" and not DROPLETQC_CONFIG.get("bam_file"):
    raise ValueError(
        "QC_RNA.DropletQC.bam_file must be set when "
        "QC_RNA.background_correction is 'soupx'."
    )


def get_dropletqc_bam(wildcards):
    bam = DROPLETQC_CONFIG.get("bam_file")
    if not bam:
        raise ValueError(
            "QC_RNA.DropletQC.bam_file must be set before running the DropletQC rule."
        )
    return bam

rule CellRangerMkref:
    output:
        fasta = os.path.join(
            config['QC_RNA']['CellRangerMkref']['output_dir'],
            config['QC_RNA']['CellRangerMkref']['genome'],
            "fasta/genome.fa"
            )
    params:
        cellranger = config['User']['cellranger'],
        outdir = config['QC_RNA']['CellRangerMkref']['output_dir'],
        genome = config['QC_RNA']['CellRangerMkref']['genome'],
        fasta = config['QC_RNA']['CellRangerMkref']['fasta'],
        genes = config['QC_RNA']['CellRangerMkref']['genes'],
    log:
        "logs/CellRangerMkref/mkref.log"
    benchmark:
        "benchmark/CellRangerMkref/mkref.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        cd {params.outdir}
        rm -r {params.genome}
        {params.cellranger} mkref \
            --genome={params.genome} \
            --fasta={params.fasta} \
            --genes={params.genes}
        """

rule CellRangerCount:
    input:
        fastq1 = f"{INPUT}/{{sample}}_S1_L001_R1_001.fastq.gz",
        fastq2 = f"{INPUT}/{{sample}}_S1_L001_R2_001.fastq.gz",
        transcriptome = os.path.join(
            config['QC_RNA']['CellRangerMkref']['output_dir'],
            config['QC_RNA']['CellRangerMkref']['genome'],
            "fasta/genome.fa"
            )
    output:
        h5 = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        pipestance = directory("QC/RNA/CellRangerCount/{sample}"),
    params:
        cellranger = config['User']['cellranger'],
        transcriptome = os.path.join(
            config['QC_RNA']['CellRangerMkref']['output_dir'],
            config['QC_RNA']['CellRangerMkref']['genome']),
        fastq_dir = INPUT,
        localmem = 64,
    threads: 8
    log:
        "logs/CellRangerCount/{sample}.log"
    benchmark:
        "benchmark/CellRangerCount/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {output.pipestance}
        {params.cellranger}  count \
            --id={wildcards.sample} \
            --output-dir QC/RNA/CellRangerCount/{wildcards.sample} \
            --transcriptome={params.transcriptome} \
            --create-bam true \
            --fastqs={params.fastq_dir} \
            --sample={wildcards.sample} \
            --localcores={threads} \
            --localmem={params.localmem}
        """

ruleorder: PlottingAnnotationsManual > PlottingAnnotationsAutomatic

rule CellbenderRemoveBackgroundRNA:
    input:
        h5 = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        h5 = "QC/RNA/{sample}/Cellbender/{sample}_cellbender.h5",
        h5_filtered = "QC/RNA/{sample}/Cellbender/{sample}_cellbender_filtered.h5",
    params:
        cuda=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['cuda'],
        epochs=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['epochs'],
        checkpoint_mins=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['checkpoint_mins'],
        projected_ambient_count_threshold=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['projected_ambient_count_threshold'],
        empty_drop_training_fraction=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['empty_drop_training_fraction'],
        expected_cells=lambda wildcards: config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'],
        total_droplets_included=lambda wildcards: 3 * config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'],
    threads: 1
    resources:
        nvidia_gpu=1
    conda:
        "../envs/cellbender.yaml"
    log:
        "logs/CellBenderRemoveBackgroundRNA/{sample}.log"
    benchmark:
        "benchmark/CellBenderRemoveBackgroundRNA/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        cellbender remove-background \
            {params.cuda} \
            --cpu-threads {threads} \
            --epochs {params.epochs} \
            --checkpoint-mins {params.checkpoint_mins} \
            --projected-ambient-count-threshold {params.projected_ambient_count_threshold} \
            --empty-drop-training-fraction {params.empty_drop_training_fraction} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets_included} \
            --input {input.h5} \
            --output {output.h5}
        rm ckpt.tar.gz 
        """

rule CellbenderToH5ad:
    input:
        cellbender_input = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        cellbender_output = "QC/RNA/{sample}/Cellbender/{sample}_cellbender_filtered.h5",
    output:
        merged = "QC/RNA/{sample}/Cellbender/{sample}_cellbender.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/CellbenderToH5ad.py",
        outdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Cellbender/"
    conda:
        "../envs/cellbender.yaml"
    log:
        "logs/CellbenderToH5ad/{sample}.log"
    benchmark:
        "benchmark/CellbenderToH5ad/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --cellbender_input {input.cellbender_input} \
        --cellbender_output {input.cellbender_output} \
        --output {output.merged}
        """

rule DropletQC:
    input:
        filtered_h5 = "QC/RNA/CellRangerCount/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        cell_qc = "QC/RNA/{sample}/DropletQC/cell_qc.tsv",
    params:
        script = f"{workflow.basedir}/scripts/QC/DropletQC.R",
        output_dir = lambda wildcards: f"QC/RNA/{wildcards.sample}/DropletQC",
        bam = get_dropletqc_bam,
        min_nf_umi = DROPLETQC_CONFIG.get("min_nf_umi", 0.6),
    threads: 4
    conda:
        "../envs/soupx_dropletqc.yaml"
    log:
        "logs/DropletQC/{sample}.log"
    benchmark:
        "benchmark/DropletQC/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        Rscript {params.script} \
            --filtered_h5 {input.filtered_h5} \
            --output_dir {params.output_dir} \
            --min_nf_umi {params.min_nf_umi} \
            --bam {params.bam}
        """

rule SoupX:
    input:
        raw_h5 = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        filtered_h5 = "QC/RNA/CellRangerCount/{sample}/outs/filtered_feature_bc_matrix.h5",
        cell_qc = "QC/RNA/{sample}/DropletQC/cell_qc.tsv",
    output:
        cell_qc = "QC/RNA/{sample}/SoupX/cell_qc.tsv",
        soup_profile = "QC/RNA/{sample}/SoupX/soup_profile.tsv",
        matrix_dir = directory("QC/RNA/{sample}/SoupX/corrected"),
    params:
        script = f"{workflow.basedir}/scripts/QC/SoupX.R",
        output_dir = lambda wildcards: f"QC/RNA/{wildcards.sample}/SoupX",
        resolution = SOUPX_CONFIG.get("resolution", 1.0),
        contaminant_genes_arg = lambda wildcards: (
            "--contaminant_genes " + " ".join(SOUPX_CONFIG.get("contaminant_genes") or [])
            if SOUPX_CONFIG.get("contaminant_genes") else ""
        ),
    threads: 4
    conda:
        "../envs/soupx_dropletqc.yaml"
    log:
        "logs/SoupX/{sample}.log"
    benchmark:
        "benchmark/SoupX/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        Rscript {params.script} \
            --raw_h5 {input.raw_h5} \
            --filtered_h5 {input.filtered_h5} \
            --cell_qc {input.cell_qc} \
            --output_dir {params.output_dir} \
            --resolution {params.resolution} \
            {params.contaminant_genes_arg}
        """

rule SoupXToH5ad:
    input:
        raw_h5 = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        soupx_dir = "QC/RNA/{sample}/SoupX/corrected",
        cell_qc = "QC/RNA/{sample}/SoupX/cell_qc.tsv",
    output:
        h5ad = "QC/RNA/{sample}/SoupX/{sample}_soupx.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/SoupXDropletQCToH5ad.py",
        soupx_dir = lambda wildcards: f"QC/RNA/{wildcards.sample}/SoupX",
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/SoupXToH5ad/{sample}.log"
    benchmark:
        "benchmark/SoupXToH5ad/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
            --raw_h5 {input.raw_h5} \
            --soupx_dir {params.soupx_dir} \
            --output {output.h5ad} \
            --sample {wildcards.sample}
        """

rule RawFilteringRNA:
    input:
        h5ad = lambda wildcards: (
            f"{INPUT}/{wildcards.sample}.h5ad" if IS_FILTERED else
            f"QC/RNA/{wildcards.sample}/SoupX/{wildcards.sample}_soupx.h5ad" if BG_CORRECTION == "soupx" else
            f"QC/RNA/{wildcards.sample}/Cellbender/{wildcards.sample}_cellbender.h5ad"
        ),
    output:
        h5ad = "QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/RawFilteringRNA.py",
        plotdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Plots",
        donor = lambda wildcards: config['samples'][wildcards.sample]['donor'],
        min_genes = config['QC_RNA']['RawFilteringRNA']['min_genes'],
        min_cells = config['QC_RNA']['RawFilteringRNA']['min_cells'],
        n_mads = config['QC_RNA']['RawFilteringRNA']['n_mads'],
        doublet_threshold = config['QC_RNA']['RawFilteringRNA']['doublet_threshold'],
        n_top_genes = config['QC_RNA']['RawFilteringRNA']['n_top_genes'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/RawFilteringRNA/{sample}.log"
    benchmark:
        "benchmark/RawFilteringRNA/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --sample {wildcards.sample} \
        --donor {params.donor} \
        --plotdir {params.plotdir} \
        --min_genes {params.min_genes} \
        --min_cells {params.min_cells} \
        --n_mads {params.n_mads} \
        {params.is_filtered} \
        --n_top_genes {params.n_top_genes} \
        --doublet_threshold {params.doublet_threshold}
        """

rule MergeSamplesAnnData:
    input:
        h5ads = expand("QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad", sample=SAMPLES.keys()),
    output:
        h5ad = "QC/RNA/Merged/MergeSamplesAnnData/merged.filtered.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/MergeSamplesAnnData.py",
        plotdir = "QC/RNA/Merged/MergeSamplesAnnData/Plots",
        samples = " ".join(SAMPLES.keys()),
        donor_key = config['General']['donor_key'],
        sample_key = config['General']['sample_key'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/MergeSamplesAnnData/merged.log"
    benchmark:
        "benchmark/MergeSamplesAnnData/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ads} \
        --output {output.h5ad} \
        --samples {params.samples} \
        --plotdir {params.plotdir} \
        --donor_key {params.donor_key} \
        --sample_key {params.sample_key} {params.is_filtered}
        """

rule BatchCorrection:
    input:
        h5ad = "QC/RNA/Merged/MergeSamplesAnnData/merged.filtered.h5ad",
        markers =  config['QC_RNA']['BatchCorrection']['markers']
    output:
        h5ad = "QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/BatchCorrection.py",
        plotdir = "QC/RNA/Merged/BatchCorrection/Plots",
        celltype_key = config['General']['celltype_key'],
        donor_key = config['General']['donor_key'],
        sample_key = config['General']['sample_key'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/BatchCorrection/merged.log"
    benchmark:
        "benchmark/BatchCorrection/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --markers {input.markers} \
        --plotdir {params.plotdir} \
        {params.is_filtered} \
        --donor_key {params.donor_key} \
        --sample_key {params.sample_key}
        """


if not IS_ANNOTATED:
    rule ManualAnnotation:
        input:
            hd5ad = "QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad",
        output:
            manual_annotation = "QC/RNA/Merged/Annotation/merged_manual_annotation.json",
        shell:
            "touch {output.manual_annotation}"


    rule CheckManualAnnotation:
        input:
            manual_annotation="QC/RNA/Merged/Annotation/merged_manual_annotation.json"
        output:
            touch("QC/RNA/Merged/Annotation/merged_manual_annotation.checked")
        run:
            # Check if the file is empty
            if os.path.getsize(input.manual_annotation) == 0:
                raise RuntimeError(f"ERROR: Annotation file {input.manual_annotation} is empty! Please manually create a JSON with 'Cluster_number':'Celltype' format")
            
            # Check if the file is valid JSON and not empty
            try:
                with open(input.manual_annotation) as f:
                    data = json.load(f)
            except json.JSONDecodeError:
                raise RuntimeError(f"ERROR: Annotation file {input.manual_annotation} is not valid JSON! Please manually create a JSON with 'Cluster_number':'Celltype' format")
            
            if not data:
                raise RuntimeError(f"ERROR: Annotation file {input.manual_annotation} contains no data!Please manually create a JSON with 'Cluster_number':'Celltype' format")
    

rule PlottingAnnotationsManual:
    input:
        h5ad = "QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad",
        manual_annotation = "QC/RNA/Merged/Annotation/merged_manual_annotation.json",
        check = "QC/RNA/Merged/Annotation/merged_manual_annotation.checked"
    output:
        h5ad = "QC/RNA/Merged/Annotation/merged.annotated.h5ad",
        csv = "QC/RNA/Merged/Annotation/merged.CellTypeAnnotations.csv",
    params:
        script = f"{workflow.basedir}/scripts/QC/PlottingAnnotations.py",
        plotdir = "QC/RNA/Merged/Annotation/Plots",
        doublets = config['QC_RNA']['PlottingAnnotations']['doublets'],
        celltype_key = config['General']['celltype_key'],
        leiden_res = config['QC_RNA']['PlottingAnnotations']['leiden_res'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/merged.log"
    benchmark:
        "benchmark/PlottingAnnotations/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --annotation_input {input.manual_annotation} \
        --annotation_output {output.csv} \
        --plotdir {params.plotdir} \
        --celltype_key {params.celltype_key} \
        --doublets {params.doublets} \
        --leiden_res {params.leiden_res} \
        --mode manual
        """

rule PlottingAnnotationsAutomatic:
    input:
        h5ad = "QC/RNA/Merged/merged.batch_corrected.h5ad",
        manual_annotation = "QC/RNA/Merged/Annotation/merged_manual_annotation.json",
    output:
        h5ad = "QC/RNA/Merged/Annotation/merged.annotated.h5ad",
        csv = "QC/RNA/Merged/Annotation/merged.CellTypeAnnotations.csv",
    params:
        script = f"{workflow.basedir}/scripts/QC/PlottingAnnotations.py",
        plotdir = "QC/RNA/Merged/Annotation/Plots",
        doublets = config['QC_RNA']['PlottingAnnotations']['doublets'],
        celltype_key = config['General']['celltype_key'],
        leiden_res = config['QC_RNA']['PlottingAnnotations']['leiden_res'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/merged.log"
    benchmark:
        "benchmark/PlottingAnnotations/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --annotation_output {output.csv} \
        --celltype_key {params.celltype_key} \
        --plotdir {params.plotdir} \
        --doublets {params.doublets} \
        --leiden_res {params.leiden_res} \
        --mode auto
        """     

rule LabelTransfer:
    input:
        h5ad = "QC/RNA/Merged/Annotation/merged.annotated.h5ad",
    output:
        h5ad = "QC/RNA/Merged/LabelTransfer/merged.label_transferred.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/LabelTransfer.py",
        plotdir = "QC/RNA/Merged/LabelTransfer/Plots",
        reference = config['QC_RNA']['LabelTransfer']['reference_h5ad'],
        reference_columns = " ".join(config['QC_RNA']['LabelTransfer']['reference_columns']),
        celltype_key = config['General']['celltype_key'],
        n_neighbors = config['QC_RNA']['LabelTransfer'].get('n_neighbors', 15),
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/LabelTransfer/merged.log"
    benchmark:
        "benchmark/LabelTransfer/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --reference {params.reference} \
        --reference_columns {params.reference_columns} \
        --celltype_key {params.celltype_key} \
        --n_neighbors {params.n_neighbors} \
        --plotdir {params.plotdir}
        """

rule TrajectoryAnalysis:
    input:
        h5ad = "QC/RNA/Merged/Annotation/merged.annotated.h5ad",
    output:
        h5ad = "QC/RNA/Merged/TrajectoryAnalysis/merged.diffusion.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/TrajectoryAnalysis.py",
        plotdir = "QC/RNA/Merged/TrajectoryAnalysis/Plots",
        diffusion_component = config['QC_RNA']['TrajectoryAnalysis']['diffusion_component'],
        root_ctype = config['QC_RNA']['TrajectoryAnalysis']['root_ctype'],
        celltype_key = config['General']['celltype_key'],
        celltype_mask = config['QC_RNA']['TrajectoryAnalysis']['celltype_mask'],
        n_genes = config['QC_RNA']['TrajectoryAnalysis']['n_genes'],
        clustering_distance = config['QC_RNA']['TrajectoryAnalysis']['clustering_distance'],
        cat_order = config['QC_RNA']['TrajectoryAnalysis']['cat_order'],
        branching = config['QC_RNA']['TrajectoryAnalysis']['branching'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/TrajectoryAnalysis/merged.log"
    benchmark:
        "benchmark/TrajectoryAnalysis/merged.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --diffusion_component {params.diffusion_component} \
        --root_ctype "{params.root_ctype}" \
        --celltype_key {params.celltype_key} \
        --celltype_mask {params.celltype_mask} \
        --n_genes {params.n_genes} \
        --clustering_distance {params.clustering_distance} \
        --cat_order {params.cat_order} \
        --plotdir {params.plotdir} \
        --branching {params.branching}
        """   
