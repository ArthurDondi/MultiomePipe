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
        expand("QC/RNA/{{sample}}/SplittingInput/split_{i}_{{sample}}.raw_feature_bc_matrix.h5", i=SPLITS),
        pngout = "QC/RNA/{sample}/Plots/{sample}.ElbowPlot.png"
    params:
        script=f"{workflow.basedir}/scripts/QC/SplittingInput.py",
        outdir=lambda wildcards: f"QC/RNA/{wildcards.sample}/SplittingInput/",
        n_splits = config['QC_RNA']['CellbenderRemoveBackgroundRNA']['n_splits']
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/SplittingInput/{sample}.log"
    benchmark:
        "benchmark/SplittingInput/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --h5in {input} \
        --pngout {output.pngout} \
        --outdir {params.outdir} \
        --sample {wildcards.sample} \
        --n_splits {params.n_splits}
        """    

rule CellbenderRemoveBackgroundRNA:
    input:
        h5 = "QC/RNA/{sample}/SplittingInput/split_{i}_{sample}.raw_feature_bc_matrix.h5"
    output:
        h5 = "QC/RNA/{sample}/Cellbender/split_{i}_{sample}_cellbender.h5",
        h5_filtered = "QC/RNA/{sample}/Cellbender/split_{i}_{sample}_cellbender_filtered.h5",
    params:
        cuda=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['cuda'],
        epochs=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['epochs'],
        checkpoint_mins=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['checkpoint_mins'],
        projected_ambient_count_threshold=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['projected_ambient_count_threshold'],
        empty_drop_training_fraction=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['empty_drop_training_fraction'],
        expected_cells=lambda wildcards: config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'],
        total_droplets_included=lambda wildcards: 3 * config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'],
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
            --empty-drop-training-fraction {params.empty_drop_training_fraction} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets_included} \
            --input {input.h5} \
            --output {output.h5}
        rm ckpt.tar.gz 
        """

rule MergingCellbenderOutput:
    input:
        cellbender_input = expand("QC/RNA/{{sample}}/SplittingInput/split_{i}_{{sample}}.raw_feature_bc_matrix.h5", i=SPLITS),
        cellbender_output = expand("QC/RNA/{{sample}}/Cellbender/split_{i}_{{sample}}_cellbender_filtered.h5", i=SPLITS),
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
        --cellbender_input {input.cellbender_input} \
        --cellbender_output {input.cellbender_output} \
        --output {output.merged}
        """

rule RawFilteringRNA:
    input:
        h5ad = lambda wildcards: f"{INPUT}/{{sample}}.h5ad" if IS_FILTERED else "QC/RNA/{sample}/Cellbender/merged_{sample}_cellbender.h5ad",
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
        donor_key = config['QC_RNA']['MergeSamplesAnnData']['donor_key'],
        sample_key = config['QC_RNA']['MergeSamplesAnnData']['sample_key'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/MergeSamplesAnnData/merged.log"
    benchmark:
        "benchmark/MergeSamplesAnnData/merged.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5ads} \
        --output {output.h5ad} \
        --samples {params.samples} \
        --plotdir {params.plotdir} \
        --donor_key {params.donor_key} \
        {params.is_filtered} \
        --sample_key {params.sample_key}
        """

rule BatchCorrection:
    input:
        h5ad = "QC/RNA/Merged/MergeSamplesAnnData/merged.filtered.h5ad",
        markers =  f"{INPUT}/marker_genes.json",
    output:
        h5ad = "QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/BatchCorrection.py",
        plotdir = "QC/RNA/Merged/BatchCorrection/Plots",
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
        donor_key = config['QC_RNA']['MergeSamplesAnnData']['donor_key'],
        sample_key = config['QC_RNA']['MergeSamplesAnnData']['sample_key'],
        is_filtered = "--is_filtered" if IS_FILTERED else ""
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
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
        leiden_res = config['QC_RNA']['PlottingAnnotations']['leiden_res'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/merged.log"
    benchmark:
        "benchmark/PlottingAnnotations/merged.benchmark.txt"
    shell:
        r"""
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
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
        leiden_res = config['QC_RNA']['PlottingAnnotations']['leiden_res'],
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/PlottingAnnotations/merged.log"
    benchmark:
        "benchmark/PlottingAnnotations/merged.benchmark.txt"
    shell:
        r"""
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
        celltype_key = config['QC_RNA']['PlottingAnnotations']['celltype_key'],
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