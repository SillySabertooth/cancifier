#since it could be hard to unite 5000 series one by one when each seria is 60000 lenght
#i divided them on 5 bathces and that process them

import pandas as pd
import numpy as np


#create an empty list of transcripts above
for i in ["aa","ab","ac","ad","ae"]:
    df = pd.read_csv("data/TCGA-4C-A93U-01A.htseq.counts.gz", compression='gzip', 
                 names = ["transcripts", "counts"], sep='\t', error_bad_lines=False).drop("counts", axis=1)
    with open("chunks."+i) as file:
        for line in file:
            line=line.rstrip()
            df_name=line.replace(".htseq.counts.gz","")
            df_new = pd.read_csv("data/"+line, compression='gzip', names = ["transcripts", df_name], sep='\t', error_bad_lines=False)
            #df_new.info()
            df = pd.merge(df, df_new)
            #print(df_new.head())
    print(df.shape)
    df.to_csv('chunks.'+i+'.tsv',sep='\t',index=False)

df1 = pd.read_csv("chunks.aa.tsv", sep="\t")
df2 = pd.read_csv("chunks.ab.tsv", sep="\t")
df3 = pd.read_csv("chunks.ac.tsv", sep="\t")
df4 = pd.read_csv("chunks.ad.tsv", sep="\t")
df5 = pd.read_csv("chunks.ae.tsv", sep="\t")
df = pd.merge(df1, df2)
df = pd.merge(df, df3)
df = pd.merge(df, df4)
df = pd.merge(df, df5)
print(df.shape)

df.to_csv('TCGA_data.tsv',sep='\t',index=False)
