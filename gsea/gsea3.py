import argparse,os,re
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join
from gseapy.scipalette import SciPalette
NbDr = SciPalette().create_colormap()
import gseapy as gp 
from gseapy import Msigdb

cluster_name='leiden_sub_1.40'
metafile='le.csv'

def readtsv(file):
    with open(file) as f:
        lst = f.read().split('\n')
    lst = [l.split('\t') for l in lst if l]
    return lst

def rename(name):
    name=re.sub(' ','_',name)
    name=re.sub(',','_',name)
    name=re.sub('/','_',name)
    return name

def main():
    #data = ad.read_h5ad('all5.h5')
    meta = pd.read_csv(metafile,index_col=0)
    #meta = pd.read_csv(metafile,sep='\t',index_col=0)
    #meta = meta[~meta[cluster_name].str.contains('|',regex=False)]
    bg = readtsv('bg.txt')
    bg = [b[0] for b in bg if b]
    print(meta)
    meta = meta[meta[cluster_name] != 'Unassigned']
    data = ad.read_h5ad('macro.h5ad')
    data = data[data.obs.index.isin(meta.index)]
    data.obs[cluster_name] = meta[cluster_name]

    print(data.obs.loc[data.obs[cluster_name] == f'SPP1 high C1QA high recruited macrophage',cluster_name])
    data.obs.loc[data.obs[cluster_name] == f'SPP1 high C1QA high recruited macrophage',cluster_name] = f'SPP1 high recruited macrophage'

    #exit()
    obs = data.obs
    #metatot = pd.read_csv('all.tsv',sep='\t',index_col=0)
    #obs = pd.concat([meta[[cluster_name]],metatot[['sample_name']]],join='inner',axis=1)
    print(obs)
    print(obs.index)
    #plotpro(pro,'Monocyte,DC,Macrophage','Monocyte_DC_Macrophage.png')
    #dic = {cluster_name:'celltype','macrotypes':'macrotype'}
    dic = {cluster_name:'celltype'}
    #data = ad.read_h5ad('all5.h5')
    #data = data[obs.index]
    #data.obs = obs
    
    data = data[        
        obs.sample_name.str.startswith('d10') |
        obs.sample_name.str.startswith('d14') |
        obs.sample_name.str.startswith('d21')
    ]

    #data.obs['sample_group'] = obs.sample_name.copy()
    #for day in ['d10','d14','d21']:
    #    for med in ['BI','Cont']:
    #        data.obs.loc[data.obs.sample_group.str.startswith(f'{day}_{med}'),'sample_group'] = f'{day}_{med}'
    
    msig = Msigdb()
    gmts = {
        'mh':msig.get_gmt(category='mh.all',dbver='2024.1.Mm'),
        #'m2':msig.get_gmt(category='m2.all',dbver='2024.1.Mm'),
        #'m3':msig.get_gmt(category='m3.all',dbver='2024.1.Mm'),
        #'m5':msig.get_gmt(category='m5.all',dbver='2024.1.Mm'),
        #'m8':msig.get_gmt(category='m8.all',dbver='2024.1.Mm')
    }

    target = [
        'HALLMARK_APOPTOSIS',
        'HALLMARK_P53_PATHWAY',
        'HALLMARK_GLYCOLYSIS',
        #'HALLMARK_INFLAMMATORY_RESPONSE',
        'HALLMARK_MTORC1_SIGNALING',
        'HALLMARK_COMPLEMENT',
        'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
        'HALLMARK_DNA_REPAIR',
        'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
    ]

    for cname in [cluster_name]:
        odir = dic[cname]
        os.makedirs(odir,exist_ok=True)
        names = data.obs[cname].unique()
        for key,gmt in gmts.items():
            oodir = os.path.join(odir,key)
            os.makedirs(oodir,exist_ok=True)
            for name in names:
                da = data[data.obs[cname] == name]
                groups=da.obs.sample_group.unique()
                sc.tl.rank_genes_groups(da,'sample_group',method='wilcoxon')
                granks = sc.get.rank_genes_groups_df(da,None)
                for day in ['d10','d14','d21']:
                        if f'{day}_BI' in groups and f'{day}_Cont' in groups:
                            for grp in ['BI','Cont']:
                                print(granks)
                                grank = granks[granks.group==f'{day}_{grp}'][['names','logfoldchanges']]
                                #sc.tl.rank_genes_groups(da,groupby='sample_group',groups=[f'{day}_BI'],reference=f'{day}_Cont')
                                #f = sc.get.rank_genes_groups_df(da,group=None)
                                #print(df)
                                #print(gmt)
                                #grank = df[['names','logfoldchanges']]
                                grank.sort_values(by=['logfoldchanges'],inplace=True,ascending=False)
                                #grank.names = grank.names.str.upper()
                                grank = grank.set_index('names')
                                print(grank)
                                res = gp.prerank(
                                    rnk = grank,
                                    gene_sets = {t:gmt[t] for t in target},
                                    threads = 12,
                                    min_size = 5,
                                    max_size = 1000,
                                    permutation_num = 10000,
                                    oudir = None,
                                    seed = 6,
                                    verbose = True,
                                )
                                terms = res.res2d.Term[1:11]
                                res.plot(terms=terms,show_ranking=True,ofname=join(oodir,f'{rename(name)}_{day}_{grp}.png'))

if __name__ == "__main__" :
    main()
