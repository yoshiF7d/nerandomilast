import anndata as ad
import pandas as pd

data = ad.read_h5ad('all5.h5')
meta = pd.read_csv('250623meta.tsv',sep='\t',index_col=0,dtype=object)
meta = meta[meta.columns[~meta.columns.str.endswith('(1)') & ~meta.columns.str.endswith('(2)')]]
data = data[data.obs.index.isin(meta.index)]
print(data.obs)
print(meta)
data.obs = meta
#data = data[data.obs['Epithelium sub_resolution0.2']!='Unassigned']
print(data)
data.write('epi2.h5ad')
