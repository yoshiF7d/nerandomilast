import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

sc.settings.dpi = 300

def rename(name):
    name=re.sub(' ','_',name)
    name=re.sub(',','_',name)
    name=re.sub('/','_',name)
    return name

filein='all5.h5'
#cluster="Epithelium sub_resolution0.2"
#metafile='250623meta.tsv'

metafile='250828meta.tsv'
cluster='leiden_res_0.20'
#rem=['16','17','18','19']
#metafile2='leiden_sub_1.40_monomacDC_metadata.tsv'

data = ad.read_h5ad(filein)
meta = pd.read_csv(metafile,sep='\t',index_col=0)
data.obs[cluster] = meta[cluster]
data.obs['sample_group'] = meta['sample_group']
#meta2 = pd.read_csv('sample_group.csv',index_col=0)
#data.obs['sample_group'] = meta2['sample_group']

#data = data[~data.obs[cluster].isin(rem)]
data = data[data.obs.sample_group!='d21_Nint']

#print(data.obs[cluster].str.startswith('Classical'))
#exit()
data = data[~data.obs[cluster].str.startswith('Classical')]
data.obs.loc[data.obs[cluster].str.endswith('monocyte'),cluster] = 'Monocyte'
data = data[data.obs[cluster] != 'Unassigned']
data = data[data.obs.sample_group != 'Unassigned']
for g in data.obs.sample_group.unique():
    print(g)
    da = data[data.obs.sample_group == g]
    print(da.obs[cluster].value_counts())
    fig,ax = plt.subplots()
    keys = np.sort(da.obs[cluster].unique())
    colors = da.uns[f'{cluster}_colors']
    cdic = dict(zip(keys,colors))
    custom = [Line2D([],[],marker='.',color=c,linestyle='None') for c in colors]
    clst = list(map(cdic.get,da.obs[cluster]))
    ax.scatter(da.obsm['X_umap'][:,0],da.obsm['X_umap'][:,1],s=0.5,lw=0,c=clst,alpha=np.clip(da.X.sum(axis=1)/2000,0,1))
    ax.legend(handles=custom,labels=list(keys),fontsize='small',bbox_to_anchor=(1.05,0.5),loc='center left')
    ax.set_title(g)
    ax.set_xlim([-7.4,19.6])
    ax.set_ylim([-10,23.3])
    ax.axis('off')
    fig.savefig(rename(g)+'.png',bbox_inches='tight',dpi=300)
    #da = data[data.obs.sample_group == g]
    #sc.pl.umap(da,color=cluster,size=4,show=False,save='_'+g+'_umap.png')
#sc.pl.umap(data,color='sample_group',size=4,show=False,save='_umap_sample.png')

#sc.pl.umap(data,color=[cluster,'sample_name'],size=2,show=False,save='_umap.png')

#data.obs = pd.merge(data.obs,meta,how='left',left_index=True,right_index=True)

