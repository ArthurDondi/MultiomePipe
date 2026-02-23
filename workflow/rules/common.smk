import pandas as pd

# Run
QC_RNA=config['Run']['QC_RNA']
QC_ATAC=config['Run']['QC_ATAC']
TRAJ = config['Run']['TRAJ']

INPUT = config['User']['input_dir']
SAMPLES = config["samples"]

IS_10X_REPO = config['IS_10X_REPO']
IS_LOCAL = config.get('IS_LOCAL', False)
IS_FILTERED = config['IS_FILTERED']
IS_ANNOTATED = config['IS_ANNOTATED']
IS_BATCHED = config['IS_BATCHED']


def get_CTYPES(samples):
    CTYPES = {}
    for sample in samples:
        ctype_file = f'{INPUT}/{sample}/cell_data.tsv'
        df_ctype = pd.read_csv(ctype_file, sep = '\t')
        ctype = list(df_ctype.barcodes.unique())
        CTYPES[sample] = ctype
    return CTYPES

def _get_macs_outputs(wildcards):
    # produce list of expected MACS narrowPeak files (one per celltype)
    celltypes = _get_celltypes_for_sample(wildcards)
    return expand("QC/ATAC/MACS/{{sample}}/{{sample}}_{celltype}_peaks.narrowPeak", celltype=CTYPES[wildcards.sample])

#CTYPES = get_CTYPES()

