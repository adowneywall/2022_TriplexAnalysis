#!/usr/bin/env python3

import pandas as pd
import os
import sys

## Script summarizes all lncRNA TTS bed files and counts number of TTS-TFO pairs for each lncRNA for each gene region

DIR2 = sys.argv[1]

import pandas as pd
import os

df_c = pd.DataFrame()
for root, dirs, files in os.walk(DIR2):
    for file in files:
        if file.endswith("intersect.bed"):
            print(file)
            df = pd.read_table(os.path.join(root,file),header=None,low_memory=False)

            gene_count = df.groupby(3).count()[1]
            gene_count = pd.DataFrame({'gene' : gene_count.index,
                                      'count' : gene_count.values})

            file_name = os.path.splitext(file)[0]
            fn = file_name.split("_")[0] + file_name.split("_")[3] + "_" + file_name.split("_")[5]

            gene_count = gene_count.rename(columns={'count' : fn})

            if not df_c.empty:
                df_c = df_c.merge(gene_count, on='gene', how='outer')
            else:
                df_c=gene_count

df_c.to_csv("TTS_CountSummary.csv")