### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets

import os
import json

ruleorder: PlottingAnnotationsManual > PlottingAnnotationsAutomatic

rule SplittingInput:
    input:
        f"{INPUT}/{{sample}}.raw_feature_bc_matrix.h5"
    output:
        expand("QC/RNA/{{sample}}/SplittingInput/split_{i}.raw_feature_bc_matrix.h5", i=SPLITS),
    params:
        script=f"{workflow.basedir}/scripts/QC/SplittingInput.py",
        outdir=lambda wildcards: f"QC/RNA/{wildcards.sample}/SplittingInput/",
        plotdir=lambda wildcards: f"QC/RNA/{wildcards.sample}/Plots",
        n_splits = config['QC_RNA']['CellbenderRemoveBackgroundRNA']['n_splits']
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/SplittingInput/{sample}.log"
    benchmark:
        "benchmark/SplittingInput/{sample}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --h5in {input} \
        --plotdir {params.plotdir} \
        --outdir {params.outdir} \
        --sample {wildcards.sample} \
        --n_splits {params.n_splits}
        """    

rule CellbenderRemoveBackgroundRNA:
    input:
        h5 = "QC/RNA/{sample}/SplittingInput/split_{i}.raw_feature_bc_matrix.h5"
    output:
        h5 = "QC/RNA/{sample}/Cellbender/split_{i}_{sample}_cellbender.h5",
    params:
        cuda=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['cuda'],
        epochs=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['epochs'],
        checkpoint_mins=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['checkpoint_mins'],
        projected_ambient_count_threshold=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['projected_ambient_count_threshold'],
        #empty_drop_training_fraction=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['empty_drop_training_fraction'],
        #expected_cells=lambda wildcards: config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'][wildcards.sample],
        #total_droplets_included=lambda wildcards: 3 * config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'][wildcards.sample],
    threads: 1
    conda:
        "../envs/cellbender.yaml"
    log:
        "logs/CellBenderRemoveBackgroundRNA/{sample}.{i}.log"
    benchmark:
        "benchmark/CellBenderRemoveBackgroundRNA/{sample}.{i}.benchmark.txt"
    shell:
        r"""
        cellbender remove-background \
            {params.cuda} \
            --cpu-threads {threads} \
            --epochs {params.epochs} \
            --checkpoint-mins {params.checkpoint_mins} \
            --projected-ambient-count-threshold {params.projected_ambient_count_threshold} \
            --input {input.h5} \
            --output {output.h5}
        rm ckpt.tar.gz 
        """

rule MergingCellbenderOutput:
    input:
        expand("QC/RNA/{{sample}}/Cellbender/split_{i}_{{sample}}_cellbender.h5", i=SPLITS),
    output:
        merged = "QC/RNA/{sample}/Cellbender/merged_{sample}_cellbender.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/MergingCellbenderOutput.py",
        outdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Cellbender/"
    conda:
        "../envs/cellbender.yaml"
    log:
        "logs/MergingCellbenderOutput/{sample}.log"
    benchmark:
        "benchmark/MergingCellbenderOutput/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input} \
        --output {output.merged}
        """

rule RawFilteringRNA:
    input:
        h5ad = lambda wildcards: f"{INPUT}/{{sample}}.h5ad" if IS_FILTERED else "QC/RNA/{sample}/Cellbender/merged_{sample}_cellbender.h5ad",
        markers =  f"{INPUT}/marker_genes.json"
    output:
        h5ad = "QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/RawFilteringRNA.py",
        plotdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Plots",
        min_genes = config['QC_RNA']['RawFilteringRNA']['min_genes'],
        min_cells = config['QC_RNA']['RawFilteringRNA']['min_cells'],
        doublet_threshold = config['QC_RNA']['RawFilteringRNA']['doublet_threshold'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/RawFilteringRNA/{sample}.log"
    benchmark:
        "benchmark/RawFilteringRNA/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --markers {input.markers} \
        --sample {wildcards.sample} \
        --plotdir {params.plotdir} \
        --min_genes {params.min_genes} \
        --min_cells {params.min_cells} \
        {params.is_filtered} \
        --doublet_threshold {params.doublet_threshold}
        """

if not IS_ANNOTATED:
    rule ManualAnnotation:
        input:
            hd5ad = "QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad",
        output:
            manual_annotation = "QC/RNA/{sample}/Annotation/{sample}_manual_annotation.json",
        shell:
            "touch {output.manual_annotation}"


    rule CheckManualAnnotation:
        input:
            manual_annotation="QC/RNA/{sample}/Annotation/{sample}_manual_annotation.json"
        output:
            touch("QC/RNA/{sample}/Annotation/{sample}_manual_annotation.checked")
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
        h5ad = "QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad",
        manual_annotation = "QC/RNA/{sample}/Annotation/{sample}_manual_annotation.json",
        check = "QC/RNA/{sample}/Annotation/{sample}_manual_annotation.checked"
    output:
        h5ad = "QC/RNA/{sample}/Annotation/{sample}_annotated.h5ad",
        csv = "QC/RNA/{sample}/Annotation/{sample}_CellTypeAnnotations.csv",
    params:
        script = f"{workflow.basedir}/scripts/QC/PlottingAnnotations.py",
        plotdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Plots",
        doublets = config['QC_RNA']['PlottingAnnotations']['doublets'],
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/{sample}.log"
    benchmark:
        "benchmark/PlottingAnnotations/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --annotation_input {input.manual_annotation} \
        --annotation_output {output.csv} \
        --sample {wildcards.sample} \
        --plotdir {params.plotdir} \
        --celltype_key {params.celltype_key} \
        --doublets {params.doublets} \
        --mode manual
        """

rule PlottingAnnotationsAutomatic:
    input:
        h5ad = "QC/RNA/{sample}/Filtering/{sample}_filtered.h5ad",
    output:
        h5ad = "QC/RNA/{sample}/Annotation/{sample}_annotated.h5ad",
        csv = "QC/RNA/{sample}/Annotation/{sample}_CellTypeAnnotations.csv",
    params:
        script = f"{workflow.basedir}/scripts/QC/PlottingAnnotations.py",
        plotdir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Plots",
        doublets = config['QC_RNA']['PlottingAnnotations']['doublets'],
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/{sample}.log"
    benchmark:
        "benchmark/PlottingAnnotations/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --annotation_output {output.csv} \
        --celltype_key {params.celltype_key} \
        --sample {wildcards.sample} \
        --plotdir {params.plotdir} \
        --doublets {params.doublets} \
        --mode auto
        """     
rule ConcatanateAnnDataObjects:
    input:
        h5ads = expand("QC/RNA/{sample}/Annotation/{sample}_annotated.h5ad", sample=SAMPLES.keys()),
        markers =  f"{INPUT}/marker_genes.json"
    output:
        h5ad = "QC/RNA/Merged/merged.annotated.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/ConcatanateAnnDataObjects.py",
        plotdir = lambda wildcards: f"QC/RNA/Merged/Plots",
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
        batch_key = config['QC_RNA']['ConcatanateAnnDataObjects']['batch_key'],
        sample_key = config['QC_RNA']['ConcatanateAnnDataObjects']['sample_key'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/ConcatanateAnnDataObjects/merged.log"
    benchmark:
        "benchmark/ConcatanateAnnDataObjects/merged.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ads} \
        --output {output.h5ad} \
        --markers {input.markers} \
        --plotdir {params.plotdir} \
        --celltype_key {params.celltype_key} \
        --batch_key {params.batch_key} \
        --sample_key {params.sample_key} \
        {params.is_filtered}
        """

rule BatchCorrection:
    input:
        h5ad = "QC/RNA/Merged/merged.annotated.h5ad",
    output:
        h5ad = "QC/RNA/Merged/merged.batch_corrected.h5ad",
        csv = "QC/RNA/Merged/CellTypeAnnotations.csv"
    params:
        script = f"{workflow.basedir}/scripts/QC/BatchCorrection.py",
        plotdir = lambda wildcards: f"QC/RNA/Merged/Plots",
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
        batch_key = config['QC_RNA']['ConcatanateAnnDataObjects']['batch_key'],
        sample_key = config['QC_RNA']['ConcatanateAnnDataObjects']['sample_key'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/BatchCorrection/merged.log"
    benchmark:
        "benchmark/BatchCorrection/merged.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ad} \
        --output {output.h5ad} \
        --annotations {output.csv} \
        --plotdir {params.plotdir} \
        --celltype_key {params.celltype_key} \
        --batch_key {params.batch_key} \
        --sample_key {params.sample_key}
        """

