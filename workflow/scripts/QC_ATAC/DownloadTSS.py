#!/usr/bin/env python3
"""
Download TSS information from pybiomart

Usage:
python DownloadTSS.py --output outs/consensus_peak_calling/sampleA_consensus_peaks.bed \
    --sample sampleA \
    --fragments path/to/fragments.tsv.gz
"""

import argparse
import pybiomart as pbm

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--tss_dataset", default="hsapiens_gene_ensembl")
    p.add_argument("--tss_host", default="http://www.ensembl.org")

    args = p.parse_args()

    dataset = pbm.Dataset(name = args.tss_dataset, host = args.tss_host)
    
    annot = dataset.query(
        attributes=[
            'chromosome_name',
            'transcription_start_site',
            'strand',
            'external_gene_name',
            'transcript_biotype']
    )

    # Formating

    # Keeping standard chroms 
    valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].astype('str')
    filt = annot['Chromosome/scaffold name'].isin(valid_chroms)
    annot = annot[filt]

    # Changing from '1' to 'chr1' format
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
    
    annot.columns = ['# Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Strand'] = annot['Strand'].replace({-1:'-', 1:'+'})

    # Only keep protein coding genes TSS (Might want to change this)
    annot = annot[annot.Transcript_type == 'protein_coding'] 

    annot.to_csv(args.output, sep='\t', index = False)

if __name__ == "__main__":
    main()
