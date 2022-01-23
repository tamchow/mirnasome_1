#!/usr/bin/env python3

import pandas as pd
import sys


sample_miRNA_df = pd.read_table(sys.argv[1], index_col='#miRNA')
sample_miRNA_df.index.name = 'miRNA'
sample_miRNA_df.drop(columns=['precursor','total','read_count'], inplace=True)
sample_id = sys.argv[2]
sample_miRNA_df.rename(columns={col: f"{col}_{sample_id}" for col in sample_miRNA_df.columns}, inplace=True)
sample_miRNA_df = sample_miRNA_df[~sample_miRNA_df.index.duplicated()]
try:
    output_miRNA_df = pd.read_csv(sys.argv[3], index_col='miRNA')
    output_miRNA_df = output_miRNA_df.join(sample_miRNA_df, how='inner', rsuffix=f'_{sample_id}')
    output_miRNA_df.to_csv(sys.argv[3])
except Exception:
    sample_miRNA_df.to_csv(sys.argv[3])
