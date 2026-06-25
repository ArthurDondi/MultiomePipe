import re
import pandas as pd

# Run
QC_RNA=config['Run']['QC_RNA']
QC_ATAC=config['Run']['QC_ATAC']
TRAJ = config['Run']['TRAJ']
LABEL_TRANSFER = config['Run'].get('LABEL_TRANSFER', False)

INPUT = config['User']['input_dir']
SAMPLES = config["samples"]

IS_10X_REPO = config['IS_10X_REPO']
IS_LOCAL = config.get('IS_LOCAL', False)
IS_FILTERED = config['IS_FILTERED']
IS_ANNOTATED = config['IS_ANNOTATED']
IS_BATCHED = config['IS_BATCHED']


# Constrain the `sample` wildcard to the configured sample names. Without this,
# `{sample}` defaults to a greedy `.+` that matches `/`, so a request for a file
# nested under a directory() output (e.g. the CellRanger pipestance file
# QC/RNA/CellRangerCount/<sample>/outs/filtered_feature_bc_matrix.h5) gets
# matched with `sample` swallowing the sub-path. That bogus binding then breaks
# input functions that look up SAMPLES[wildcards.sample] (KeyError).
wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),


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

