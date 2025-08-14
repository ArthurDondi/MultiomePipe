# Download raw or filtered data

rule GetRawData:
    output:
        f"{INPUT}/raw/{{sample}}.raw_feature_bc_matrix.h5"
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
