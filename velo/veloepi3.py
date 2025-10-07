import scanpy as sc
import anndata as ad
import pandas as pd
import os
import scvelo as scv
from multiprocessing import Pool,set_start_method

projects = [
    #['F8484','d0','d8_1','d8_2'],
    ['F8655','d10_BI_1','d10_BI_2','d10_Cont_1','d10_Cont_2'],
    ['F8668','d14_BI_1','d14_BI_2','d14_Cont_1','d14_Cont_2'],
    ['F8492','d21_BI_1','d21_BI_2','d21_Cont_1','d21_Cont_2']
]
cluster = 'Epithelium sub_resolution0.2'
metafile = '250623meta.tsv'

class GOI:
    def __init__(self,name,proj,data):
        self.name = name
        self.proj = proj
        self.data = data
    def print(self):
        print(self.name)
        print(self.proj)
        print(self.data)

def go(goi):
    loom = sc.read(os.path.join('../../',goi.proj,goi.name,'velocyto',goi.name+'.loom'))
    merge = scv.utils.merge(goi.data,loom)
    merge.write('loom_'+goi.name+'.h5ad')

def main():
    adata = ad.read_h5ad('epi.h5ad')
    #meta = pd.read_csv(metafile,sep='\t',index_col=0)
    #obs = pd.concat([meta[[cluster]],adata.obs[['sample_name']]],join='inner',axis=1)
    #adatap = adata[obs.index]
    #ndata = ad.AnnData(X=adatap.X,obs=obs,var=adatap.var,uns={'umap':adatap.uns['umap']},obsm=adatap.obsm,varm=adatap.varm)
    
    gois = []
    for p in projects:
        for name in p[1:]:
            print(p[0]+' '+name)
            #print(ndata[adata.obs.sample_name==name])
            gois.append(GOI(name,p[0],adata[adata.obs.sample_name==name].copy()))
    with Pool(processes=12) as pool:
       print('go')
       pool.map(go,gois) 

if __name__ == "__main__" :
    set_start_method('spawn')
    main()
