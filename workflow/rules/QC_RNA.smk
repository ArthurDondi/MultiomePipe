### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets

import os
import json

# Background correction method: "soupx" (default, DropletQC + SoupX) or "cellbender"
BG_CORRECTION = config.get("QC_RNA", {}).get("background_correction", "soupx")
if BG_CORRECTION not in {"cellbender", "soupx"}:
    raise ValueError(
        f"Unsupported QC_RNA.background_correction '{BG_CORRECTION}'. "
        "Choose 'cellbender' or 'soupx'."
    )

# Batch correction method: "harmony" (default, CPU) or "scvi"/"scanvi" (GPU).
# scvi/scanvi run on the `gpu` SLURM queue and use the scvi-tools conda env.
_BC_CFG = config.get("QC_RNA", {}).get("BatchCorrection", {})
BC_METHOD = _BC_CFG.get("method", "harmony").lower()
if BC_METHOD not in {"harmony", "scvi", "scanvi"}:
    raise ValueError(
        f"Unsupported QC_RNA.BatchCorrection.method '{BC_METHOD}'. "
        "Choose 'harmony', 'scvi' or 'scanvi'."
    )
BC_USE_GPU = BC_METHOD in {"scvi", "scanvi"}
BC_ENV = "../envs/scvi.yaml" if BC_USE_GPU else "../envs/scverse.yaml"
# Wall-clock for the BatchCorrection job (minutes). GPU runs default to 4h on
# the gpu queue; harmony keeps the historical 6h. Override via config.
BC_RUNTIME = int(_BC_CFG.get("runtime", 240 if BC_USE_GPU else 360))


def _bc_scvi_args(cfg):
    """CLI flags for scVI/scANVI hyperparameters and label sources."""
    if not BC_USE_GPU:
        return ""
    scvi = cfg.get("scvi", {}) or {}
    parts = [f"--method {BC_METHOD}"]
    parts.append(f"--n_hvg {scvi.get('n_hvg', 2000)}")
    parts.append(f"--n_latent {scvi.get('n_latent', 30)}")
    if scvi.get("max_epochs") is not None:
        parts.append(f"--scvi_max_epochs {scvi['max_epochs']}")
    if BC_METHOD == "scanvi":
        scanvi = cfg.get("scanvi", {}) or {}
        parts.append(f"--scanvi_max_epochs {scanvi.get('max_epochs', 20)}")
        parts.append(f"--unlabeled_category '{scanvi.get('unlabeled_category', 'Unknown')}'")
        sources = scanvi.get("label_sources", []) or []
        if sources:
            ds = " ".join(str(s["dataset"]) for s in sources)
            refs = " ".join(str(s["reference_h5ad"]) for s in sources)
            cols = " ".join(str(s["label_column"]) for s in sources)
            parts.append(f"--label_datasets {ds}")
            parts.append(f"--label_refs {refs}")
            parts.append(f"--label_columns {cols}")
    return " ".join(parts)

rule CellRangerMkref:
    input:
        # Track the source FASTA + GTF as inputs (not params) so that editing
        # either one — e.g. swapping in a GTF that actually annotates chrM —
        # invalidates the built reference and forces a rebuild. As params these
        # were invisible to Snakemake, so a changed GTF left the existing
        # cellranger_ref untouched and CellRangerCount silently kept counting
        # against a stale reference (missing chrM/MT genes in features.tsv even
        # though the FASTA still had the chrM contig, hence chrM reads in the BAM).
        fasta = config['QC_RNA']['CellRangerMkref']['fasta'],
        genes = config['QC_RNA']['CellRangerMkref']['genes'],
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
    threads: 8
    resources:
        # Building a STAR index for a full human primary assembly peaks well
        # above the 8 GB / 2h profile defaults; too little of either kills
        # mkref mid-build and leaves the broken mkref_<genome> pipestance dir.
        mem_mb = 64000,
        runtime = 480,   # 8h
    log:
        "logs/CellRangerMkref/mkref.log"
    benchmark:
        "benchmark/CellRangerMkref/mkref.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        # cellranger mkref writes into the current working directory and has no
        # output-dir flag, so resolve the inputs to absolute paths before cd.
        fasta=$(realpath "{input.fasta}")
        genes=$(realpath "{input.genes}")
        mkdir -p {params.outdir}
        cd {params.outdir}
        # Remove any previous reference AND a leftover mkref pipestance dir.
        # cellranger mkref builds inside a working dir named mkref_<genome>; an
        # interrupted run (timeout/OOM/cancel) leaves a partial one behind, and
        # since it is not a declared output Snakemake never cleans it, so the
        # retry aborts with "is not a pipestance directory".
        rm -rf {params.genome} mkref_{params.genome}
        {params.cellranger} mkref \
            --genome={params.genome} \
            --fasta="$fasta" \
            --genes="$genes" \
            --nthreads={threads} \
            --memgb=$(( {resources.mem_mb} / 1000 - 8 ))
        """

rule CellRangerCount:
    input:
        # FASTQs are pre-existing raw inputs living in a per-sample directory
        # (config: samples.<sample>.fastq_dir). They are not produced by the
        # workflow, so only the transcriptome is tracked as a rule dependency.
        transcriptome = os.path.join(
            config['QC_RNA']['CellRangerMkref']['output_dir'],
            config['QC_RNA']['CellRangerMkref']['genome'],
            "fasta/genome.fa"
            )
    output:
        h5 = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        # Declared explicitly because DropletQC / SoupX consume it directly;
        # otherwise Snakemake can only reach it through the pipestance
        # directory() output, which forces a greedy {sample} match.
        filtered_h5 = "QC/RNA/CellRangerCount/{sample}/outs/filtered_feature_bc_matrix.h5",
        bam = "QC/RNA/CellRangerCount/{sample}/outs/possorted_genome_bam.bam",
        pipestance = directory("QC/RNA/CellRangerCount/{sample}"),
    params:
        cellranger = config['User']['cellranger'],
        transcriptome = os.path.join(
            config['QC_RNA']['CellRangerMkref']['output_dir'],
            config['QC_RNA']['CellRangerMkref']['genome']),
        # Per-sample FASTQ directory and CellRanger --sample prefix, with the
        # global input_dir / wildcard sample as fallbacks for older configs.
        fastq_dir = lambda wildcards: SAMPLES[wildcards.sample].get('fastq_dir', INPUT),
        cr_sample = lambda wildcards: SAMPLES[wildcards.sample].get('cellranger_sample', wildcards.sample),
        localmem = 64,
    threads: 24
    resources:
        mem_mb = 72000,
        runtime = 2880,   # 24h (mediumq)
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
            --sample={params.cr_sample} \
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
        nvidia_gpu=1,
        mem_mb = 32000,
        runtime = 720,    # 12h
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
        raw_h5      = "QC/RNA/CellRangerCount/{sample}/outs/raw_feature_bc_matrix.h5",
        filtered_h5 = "QC/RNA/CellRangerCount/{sample}/outs/filtered_feature_bc_matrix.h5",
        bam         = "QC/RNA/CellRangerCount/{sample}/outs/possorted_genome_bam.bam",
    output:
        cell_qc            = "QC/RNA/{sample}/DropletQC/cell_qc.tsv",
        cellranger_compare = "QC/RNA/{sample}/DropletQC/cellranger_comparison.tsv",
        barcode_compare    = "QC/RNA/{sample}/DropletQC/barcode_comparison.tsv",
    params:
        script             = f"{workflow.basedir}/scripts/QC/DropletQC.R",
        output_dir         = lambda wildcards: f"QC/RNA/{wildcards.sample}/DropletQC",
        assay_type         = config['QC_RNA']['DropletQC']['assay_type'],
        seurat_resolution  = config['QC_RNA']['DropletQC'].get('seurat_resolution', 0.5),
        min_umi_for_nf     = config['QC_RNA']['DropletQC'].get('min_umi_for_nf', 100),
    threads: 4
    resources:
        mem_mb = 32000,
        runtime = 360,    # 6h
    conda:
        "../envs/dropletqc.yaml"
    log:
        "logs/DropletQC/{sample}.log"
    benchmark:
        "benchmark/DropletQC/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        Rscript {params.script} \
            --raw_h5           {input.raw_h5} \
            --filtered_h5      {input.filtered_h5} \
            --bam              {input.bam} \
            --output_dir       {params.output_dir} \
            --assay_type       {params.assay_type} \
            --seurat_resolution {params.seurat_resolution} \
            --min_umi_for_nf   {params.min_umi_for_nf} \
            --cores            {threads}
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
        resolution = config['QC_RNA']['SoupX']['resolution'],
        contaminant_genes_arg = lambda wildcards: (
            "--contaminant_genes " + " ".join(config['QC_RNA']['SoupX']['contaminant_genes'])
            if config['QC_RNA']['SoupX']['contaminant_genes'] else ""
        ),
    threads: 4
    resources:
        mem_mb = 32000,
        runtime = 360,    # 6h
    conda:
        "../envs/soupx.yaml"
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
        dataset = lambda wildcards: config['samples'][wildcards.sample].get('dataset', 'NA'),
        min_genes = config['QC_RNA']['RawFilteringRNA']['min_genes'],
        min_cells = config['QC_RNA']['RawFilteringRNA']['min_cells'],
        n_mads = config['QC_RNA']['RawFilteringRNA']['n_mads'],
        doublet_threshold = config['QC_RNA']['RawFilteringRNA']['doublet_threshold'],
        n_top_genes = config['QC_RNA']['RawFilteringRNA']['n_top_genes'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    resources:
        # Scrublet (sc.pp.scrublet) ~doubles the matrix to simulate doublets,
        # so peak RAM scales with cell count; the 8 GB profile default OOM-kills
        # large samples (e.g. cell-line organoids). 32 GB covers the big ones.
        mem_mb = 32000,
        runtime = 240,    # 4h
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/RawFilteringRNA/{sample}.log"
    benchmark:
        "benchmark/RawFilteringRNA/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        python -u -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --sample {wildcards.sample} \
        --donor {params.donor} \
        --dataset {params.dataset} \
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
        dataset_key = config['General'].get('dataset_key', 'dataset'),
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    resources:
        # Peak RAM scales with the TOTAL cell count across all samples: the script
        # reads every per-sample .h5ad into memory at once, holds them while
        # building the concatenated copy, then runs PCA/neighbors/UMAP/Leiden on
        # the full merged dataset. The 8 GB profile default OOM-kills this silently
        # (empty log). 128 GB gives ample headroom for large multi-sample runs;
        # tune down if your dataset is small. 6h fits within shortq (12h).
        mem_mb = 128000,
        runtime = 360,    # 6h
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
        --sample_key {params.sample_key} \
        --dataset_key {params.dataset_key} {params.is_filtered}
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
        dataset_key = config['General'].get('dataset_key', 'dataset'),
        is_filtered = "--is_filtered" if IS_FILTERED else "",
        # method switch (+ scvi/scanvi hyperparameters & label sources); empty
        # for harmony, which then runs with the script's default method.
        method_args = _bc_scvi_args(_BC_CFG),
    resources:
        mem_mb = 32000,
        runtime = BC_RUNTIME,
        # Route scvi/scanvi to the GPU queue; harmony stays on the default
        # (shortq) queue. gres/nvidia_gpu are falsy (ignored) when not on GPU.
        slurm_partition = "gpu" if BC_USE_GPU else "shortq",
        qos = "gpu" if BC_USE_GPU else "shortq",
        gres = "gpu:1" if BC_USE_GPU else "",
        nvidia_gpu = 1 if BC_USE_GPU else 0,
    conda:
        BC_ENV
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
        --sample_key {params.sample_key} \
        --dataset_key {params.dataset_key} \
        {params.method_args}
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
    resources:
        # Re-loads the full merged (all-cells) dataset to render annotation
        # UMAPs; the 8 GB default risks OOM on large runs. Matches BatchCorrection.
        mem_mb = 32000,
        runtime = 240,    # 4h
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
    resources:
        # Re-loads the full merged (all-cells) dataset to render annotation
        # UMAPs; the 8 GB default risks OOM on large runs. Matches BatchCorrection.
        mem_mb = 32000,
        runtime = 240,    # 4h
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
    resources:
        # Loads the merged dataset AND the reference atlas, then builds a joint
        # kNN graph for label transfer — peak RAM exceeds the merged object alone.
        mem_mb = 64000,
        runtime = 360,    # 6h
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
    resources:
        # Diffusion-map / pseudotime on the merged dataset (eigendecomposition of
        # the cell-cell graph) is memory-heavy; the 8 GB default OOMs on large runs.
        mem_mb = 64000,
        runtime = 360,    # 6h
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
