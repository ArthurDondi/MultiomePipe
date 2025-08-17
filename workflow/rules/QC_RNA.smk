### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets

#if config['QC_RNA']['raw_data']:

rule RawFilteringRNA:
    input:
        f"{INPUT}/raw/{{sample}}.raw_feature_bc_matrix.h5"
    output:
        hd5 = "QC/RNA/{sample}/{sample}_filtered.hd5mu",
        clusters = "QC/RNA/{sample}/{sample}_scanpy_clusters.tsv",
        mtx = "QC/RNA/{sample}/10x_for_SoupX/matrix.mtx.gz"
    params:
        script=f"{workflow.basedir}/scripts/QC/RawFilteringRNA.py",
        outdir="QC/RNA/",
        min_genes=config['QC_RNA']['min_genes'],
        min_cells=config['QC_RNA']['min_cells']
    conda:
        "../envs/scverse.yaml"
    log:
        "logs/RawFilteringRNA/{sample}.log"
    benchmark:
        "benchmark/RawFilteringRNA/{sample}.benchmark.txt"
    shell:
        r"""
        python -W ignore {params.script} \
        --input {input} \
        --output {output.hd5} \
        --sample {wildcards.sample} \
        --outdir {params.outdir}/{wildcards.sample} \
        --min_genes {params.min_genes} \
        --min_cells {params.min_cells} 
        """
        

rule RemoveAmbientAndDoublets:
    input:
        clusters = "QC/RNA/{sample}/{sample}_scanpy_clusters.tsv",
        filtered_mtx = "QC/RNA/{sample}/10x_for_SoupX/matrix.mtx.gz",
        raw_mtx = f"{INPUT}/raw/{{sample}}/raw_feature_bc_matrix/matrix.mtx.gz"
    output:
        rds = "QC/RNA/{sample}/{sample}_processed.rds"
    params:
        script=f"{workflow.basedir}/scripts/QC/RemoveAmbientAndDoublets.R",
        filtered_dir="QC/RNA/",
        raw_dir = f"{INPUT}/raw"
    conda:
        "../envs/RemoveAmbientAndDoublets.yaml"
    log:
        "logs/RemoveAmbientAndDoublets/{sample}.log"
    benchmark:
        "benchmark/RemoveAmbientAndDoublets/{sample}.benchmark.txt"
    shell:
        r"""
        Rscript {params.script} \
        --raw_dir {params.raw_dir}/{wildcards.sample}/raw_feature_bc_matrix/ \
        --filtered_dir {params.filtered_dir}/{wildcards.sample}/10x_for_SoupX \
        --clusters {input.clusters} \
        --output {output.rds}
        """