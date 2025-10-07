import scanpy as sc
import anndata as ad
import pandas as pd
import os,re,glob,shutil,argparse
import scvelo as scv
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

samples = {
#    'd0':['d0'],
#    'd8':['d8_1','d8_2'],
    'd10_BI':['d10_BI_1','d10_BI_2'],'d10_Cont':['d10_Cont_1','d10_Cont_2'],
    'd14_BI':['d14_BI_1','d14_BI_2'],'d14_Cont':['d14_Cont_1','d14_Cont_2'],
    'd21_BI':['d21_BI_1','d21_BI_2'],'d21_Cont':['d21_Cont_1','d21_Cont_2'],
}
cluster = 'Epithelium sub_resolution0.2'
metafile = '250623meta.tsv'

for name,v in samples.items():
    das=[]
    for i,s in enumerate(v): 
        da = ad.read_h5ad(join('result',s,f'{s}.h5ad'))
        da.obs.index += ':' + str(i)
        #da.obs = da.obs[['G2M_score','S_score']]
        das.append(da)
    data = ad.concat(das)
    fig,ax = plt.subplots()
    ax.set_xlim([-4,3])
    ax.set_ylim([-4,3])
    x = data.obsm['X_umap'][:,0]
    y = data.obsm['X_umap'][:,1]
    z = data.obs.S_score
    minz,maxz=-0.05,0.4
    znorm = (z - minz)/(maxz-minz)
    ax.scatter(x,y,c='Blue',alpha=np.clip(znorm,0,1),label='S_score')
    z = data.obs.G2M_score
    minz,maxz=-0.05,1
    znorm = (z - minz)/(maxz-minz)
    ax.scatter(x,y,c='Orange',alpha=np.clip(znorm,0,1),label='G2M_score')
    handles=[]
    handles.append(mpatches.Patch(color='Blue',label='S_score'))
    handles.append(mpatches.Patch(color='Orange',label='G2M_score'))
    ax.legend(handles=handles)
    fig.savefig(f'{name}_cycle.png')
    #scv.pl.scatter(data, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95],save=name+'_cycle.png',show=False)
