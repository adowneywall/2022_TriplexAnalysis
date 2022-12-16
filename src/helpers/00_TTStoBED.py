#!/usr/bin/env python3

import pandas as pd
import os
import sys

## Script loops through each lncRNA output TTS text and splits the chromosome 
 # coordinate file and generates a gene region centric bed file

DIR = sys.argv[1]

for root, dirs, files in os.walk(DIR):
    for file in files:
        if file.endswith('.txt'):
            file_name = os.path.splitext(os.path.basename(file))[0] + '.bed'
            print(file_name)

            df = pd.read_table(os.path.join(root,file))
            df_genome = pd.DataFrame(df['DNA-ID'].str.split(':|-',))

            chr = [item[0] for item in df_genome['DNA-ID']]
            start = [item[1] for item in df_genome['DNA-ID']]
            stop = [item[2] for item in df_genome['DNA-ID']]

            pair_length = df['RBS-end']-df['RBS-start']

            df_rv = pd.DataFrame({'chr' : chr,
                          'start' : start,
                          'stop' : stop,
                          'lncRNA' : df['# RNA-ID'],
                          'lncRNA_start' : df['RBS-start'],
                          'lncRNA_end' : df['RBS-end'],
                          'score' : df['Score'],
                          'error_rate' : df['Error-rate'],
                          'g_rate' : df['Guanine-rate'],
                          'lnc_length' : pair_length})
            df_rv.to_csv(file_name, sep ='\t',header=None,index=False)