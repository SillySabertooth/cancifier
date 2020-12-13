import pandas as pd
import numpy as np

# upload gathered data
name="TCGA_data"
df = pd.read_csv(name+".tsv", sep='\t')

#adjust data to avoid -inf when log2 transf
df = df.set_index('transcripts')
#pseudo_count=0.5
#df1 = df.replace(0, pseudo_count) #deprecated way
df1 = np.log2(df+1)
print(df1.shape)
df1.to_csv(name+'.log2.tsv',sep='\t',index=False)

#quantile normalization
rank_mean = df1.stack().groupby(df1.rank(method='first').stack().astype(int)).mean()
df2 = df1.rank(method='min').stack().astype(int).map(rank_mean).unstack()
print(df2.shape)
df2.to_csv(name+'quant.log2.tsv',sep='\t',index=False)

#swithing from Ensembl_id to GeneSymbols -> match with gene names gathered from biomart (the ensembl resource)
df_anno = df2.reset_index()
ids=pd.read_csv("mart_export_2_genes.txt", delimiter= "\t")
ids.columns = ['transcripts', 'NCBI_symbol']
ids[['tr','tmp_2']] = ids.transcripts.str.split(".", expand=True)
ids = ids.drop(["transcripts","tmp_2"], axis=1)
df_anno[['tr','tmp']] = df_anno.transcripts.str.split(".", expand=True)
df_anno = df_anno.drop(["transcripts","tmp"], axis=1)

anno=pd.merge(df_anno,ids) #,how='left').fillna(0)
anno=anno.drop_duplicates()

print(anno.shape)
anno.to_csv(name+'annotated.quant.log2.tsv',sep='\t',index=False)

#removing duplicated probes by getting the probe with max median among the same probes (median counted among the samples per probe)
anno = anno.set_index('tr')
anno['median']=anno.drop("NCBI_symbol", axis=1).apply(lambda row: np.median(row), axis=1)
anno_simple = anno.sort_values('median', ascending=False).drop_duplicates(['NCBI_symbol'])

print(anno_simple.shape)
anno_simple.to_csv(name+'no_dup.annotated.quant.log2.tsv',sep='\t',index=False)

#getting the top 15000 most expressed genes
anno_simple['mean']=anno_simple.drop(["NCBI_symbol","median"], axis=1).apply(lambda row: np.mean(row), axis=1)
anno_simple.sort_values('mean', ascending=False)
anno_top = anno_simple[:15000]
anno_top.to_csv(name+'top15k.no_dup.annotated.quant.log2.tsv',sep='\t',index=False)
print(anno_top.shape, anno_top.tail())

#weird adding of batches per sample to be able to check the batch effect in the future
trans = anno_top.drop(["median","mean"], axis=1).set_index('NCBI_symbol').T #with T - transpose!!
trans = trans.reset_index()
names = trans.columns.tolist()
names[names.index('index')] = 'sample_id'
trans.columns = names

bbatch = pd.read_csv("batches.samples_TCGA", sep="\t")

new =pd.merge(trans,bbatch)
new = new.set_index("sample_id")
new = new.T
#new
new = new.reset_index()
names = new.columns.tolist()
names[names.index('index')] = 'NCBI_gene'
new.columns = names

# final dataset for TCGA (it will be pushed to Z-transfomration and diff express gene analysis)
new.to_csv(name+'.FOR_MERGE_with_BATCH.top15k.no_dup.annotated.quant.log2.tsv', sep='\t', index = False)

