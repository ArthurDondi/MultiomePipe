### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets


rule ElbowPlot:
    input:
        f"{INPUT}/raw/{{sample}}.raw_feature_bc_matrix.h5",
    output:
        expand("QC/RNA/{{sample}}/Subsampling/split_{i}.raw_feature_bc_matrix.h5", i=Is),
        png = "QC/RNA/{sample}/Plots/{sample}_ElbowPlot.png",
        #"QC/RNA/{sample}/Subsampling/split_{i}.raw_feature_bc_matrix.h5"
    params:
        script=f"{workflow.basedir}/scripts/QC/ElbowPlot.py",
        outdir="QC/RNA",
        n_splits = config['QC_RNA']['CellbenderRemoveBackgroundRNA']['n_splits']
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/ElbowPlot/{sample}.log"
    benchmark:
        "benchmark/ElbowPlot/{sample}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --h5in {input} \
        --pngout {output.png} \
        --outdir {params.outdir}/{wildcards.sample}/Subsampling/ \
        --n_splits {params.n_splits}
        """    

rule CellbenderRemoveBackgroundRNA:
    input:
        h5 = "QC/RNA/{sample}/Subsampling/split_{i}.raw_feature_bc_matrix.h5"
    output:
        h5 = "QC/RNA/{sample}/Cellbender/split_{i}_{sample}_cellbender.h5",
    params:
        cuda=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['cuda'],
        epochs=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['epochs'],
        checkpoint_mins=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['checkpoint_mins'],
        projected_ambient_count_threshold=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['projected_ambient_count_threshold'],
        empty_drop_training_fraction=config['QC_RNA']['CellbenderRemoveBackgroundRNA']['empty_drop_training_fraction'],
        expected_cells=lambda wildcards: config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'][wildcards.sample],
        total_droplets_included=lambda wildcards: 3 * config['QC_RNA']['CellbenderRemoveBackgroundRNA']['expected-cells'][wildcards.sample],
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
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets_included} \
            --projected-ambient-count-threshold {params.projected_ambient_count_threshold} \
            --empty-drop-training-fraction {params.empty_drop_training_fraction} \
            --input {input.h5} \
            --output {output.h5} \
            > {log} 2>&1
        rm ckpt.tar.gz 
        """

rule MergingCellbenderOutput:
    input:
        expand("QC/RNA/{{sample}}/Cellbender/split_{i}_{{sample}}_cellbender.h5", i=Is),
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
        h5 = "QC/RNA/{sample}/Cellbender/merged_{sample}_cellbender.h5ad",
    output:
        hd5mu = "QC/RNA/{sample}/Filtering/{sample}_filtered.hd5ad",
    params:
        script = f"{workflow.basedir}/scripts/QC/RawFilteringRNA.py",
        dir = lambda wildcards: f"QC/RNA/{wildcards.sample}/Filtering/",
        min_genes = config['QC_RNA']['min_genes'],
        min_cells = config['QC_RNA']['min_cells']
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/RawFilteringRNA/{sample}.log"
    benchmark:
        "benchmark/RawFilteringRNA/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input.h5} \
        --output {output.hd5mu} \
        --sample {wildcards.sample} \
        --outdir {params.dir} \
        --min_genes {params.min_genes} \
        --min_cells {params.min_cells} 
        """
        