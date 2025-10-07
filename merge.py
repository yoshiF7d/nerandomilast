import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
projects = [
    ['F8484','d0','d8_1','d8_2'],
    ['F8655','d10_BI_1','d10_BI_2','d10_Cont_1','d10_Cont_2'],
    ['F8668','d14_BI_1','d14_BI_2','d14_Cont_1','d14_Cont_2'],
    ['F8492','d21_BI_1','d21_BI_2','d21_Nint_1','d21_Nint_2','d21_Cont_1','d21_Cont_2']
]

def main():
    datas={}
    names=[]
    for p in projects:
        for name in p[1:]:
            names.append(name)
            data = sc.read_10x_h5(f'../{p[0]}/{name}/outs/filtered_feature_bc_matrix.h5')
            data.var_names_make_unique()
            data.var['mt'] = data.var_names.str.startswith('mt-')
            data.var['ribo'] = data.var_names.str.startswith(('Rps','Rpl'))
            sc.pp.calculate_qc_metrics(data,qc_vars=['mt','ribo'],inplace=True,log1p=True)
            sc.pp.filter_cells(data, min_genes=100)
            sc.pp.filter_genes(data, min_cells=3)
            data = data[data.obs.pct_counts_mt<20,:]
            datas[name] = data
            sc.pp.pca(data)
    
    data = ad.concat(datas,label='sample_name',index_unique='_')
    data.obs_names_make_unique()
    #12/1 fixed
    #sce.pp.harmony_integrate(data,'sample_name')
    #data.obsm['X_pca'] = data.obsm['X_pca_harmony']
    #print(data.obs)
    fig,ax = plt.subplots(2,2)
    sc.pl.violin(data,'n_genes_by_counts',stripplot=False,show=False,ax=ax[0,0])
    ax[0,0].set_ylim([0,10000])
    sc.pl.violin(data,'total_counts',stripplot=False,show=False,ax=ax[0,1])
    ax[0,1].set_ylim([0,100000])
    sc.pl.violin(data,'pct_counts_mt',stripplot=False,show=False,ax=ax[1,0])
    ax[1,0].set_ylim([0,20])
    ax[1,0].set_ylim([0,20])
    sc.pl.scatter(data,'total_counts','n_genes_by_counts',color='pct_counts_mt',show=False,ax=ax[1,1])
    ax[1,1].set_xlim([0,100000])
    ax[1,1].set_ylim([0,10000])
    fig.savefig('violin.png',bbox_inches='tight',dpi=300)
    data.write('all.h5')

if __name__ == "__main__" :
    main()
