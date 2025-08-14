### Rules for GEX data QC
### First applies rough filters with a pause to manually check
### Then remove doublets

#if config['QC_RNA']['raw_data']:

rule RawFilteringRNA:
    input:
        f"{INPUT}/raw/{{sample}}.raw_feature_bc_matrix.h5"
    output:
        hd5 = "QC/RNA/{sample}_filtered.hd5mu"
    params:
        script=f"{workflow.basedir}/scripts/QC/RawFilteringRNA.py",
        plotdir="QC/RNA/Plots",
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
        --plotdir {params.plotdir} \
        --min_genes {params.min_genes} \
        --min_cells {params.min_cells} 
        """