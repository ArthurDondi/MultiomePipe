# Download raw or filtered data

rule GetRawData:
    output:
        f"{INPUT}/{{sample}}.raw_feature_bc_matrix.h5"
    params:
        url=lambda wildcards: SAMPLES[wildcards.sample]
    log:
        "logs/GetRawData/{sample}.log"
    benchmark:
        "benchmark/GetRawData/{sample}.benchmark.txt"
    shell:
        r"""
        wget -O {output} {params.url}_raw_feature_bc_matrix.h5
        """

rule GetRawMatrix:
    output:
        f"{INPUT}/{{sample}}/raw_feature_bc_matrix/matrix.mtx.gz",
        f"{INPUT}/{{sample}}/raw_feature_bc_matrix/barcodes.tsv.gz",
        f"{INPUT}/{{sample}}/raw_feature_bc_matrix/features.tsv.gz",
    params:
        url=lambda wildcards: SAMPLES[wildcards.sample],
        outdir=f"{INPUT}/"
    log:
        "logs/GetRawMatrix/{sample}.log"
    benchmark:
        "benchmark/GetRawMatrix/{sample}.benchmark.txt"
    shell:
        r"""
        wget -O {params.outdir}/{wildcards.sample}_raw_feature_bc_matrix.tar.gz {params.url}_raw_feature_bc_matrix.tar.gz
        tar -xvzf {params.outdir}/{wildcards.sample}_raw_feature_bc_matrix.tar.gz -C {params.outdir}/{wildcards.sample}/
        """

rule GetFilteredData:
    output:
        f"{INPUT}/filtered/{{sample}}.filtered_feature_bc_matrix.h5"
    params:
        url=lambda wildcards: SAMPLES[wildcards.sample]
    log:
        "logs/GetFilteredData/{sample}.log"
    benchmark:
        "benchmark/GetFilteredData/{sample}.benchmark.txt"
    shell:
        r"""
        wget -O {output} {params.url}_filtered_feature_bc_matrix.h5
        """

rule GetATACData:
    output:
        peak_annotation = f"{INPUT}/ATAC/{{sample}}.atac_peak_annotation.tsv",
        fragments = f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz",
        fragments_table = f"{INPUT}/ATAC/{{sample}}.atac_fragments.tsv.gz.tbi"
    params:
        url=lambda wildcards: SAMPLES[wildcards.sample]
    log:
        "logs/GetATACData/{sample}.log"
    benchmark:
        "benchmark/GetATACData/{sample}.benchmark.txt"
    shell:
        r"""
        wget -O {output.peak_annotation} {params.url}.atac_peak_annotation.tsv
        wget -O {output.fragments} {params.url}.atac_fragments.tsv.gz
        wget -O {output.fragments_table} {params.url}.atac_fragments.tsv.gz.tbi
        """
