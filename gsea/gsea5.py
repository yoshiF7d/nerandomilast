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
        #'mh':msig.get_gmt(category='mh.all',dbver='2024.1.Mm'),
        #'m2':msig.get_gmt(category='m2.all',dbver='2024.1.Mm'),
        #'m3':msig.get_gmt(category='m3.all',dbver='2024.1.Mm'),
        'm5':msig.get_gmt(category='m5.all',dbver='2025.1.Mm'),
        #'m8':msig.get_gmt(category='m8.all',dbver='2024.1.Mm')
    }

    target = [
        'GOBP_CELLULAR_RESPONSE_TO_CAMP',
        'GOBP_REGULATION_OF_CAMP_DEPENDENT_PROTEIN_KINASE_ACTIVITY',
        'GOBP_REGULATION_OF_CAMP_MEDIATED_SIGNALING',
        'GOBP_REGULATION_OF_CAMP_PKA_SIGNAL_TRANSDUCTION',
        'GOBP_POSITIVE_REGULATION_OF_CAMP_MEDIATED_SIGNALING',
        'GOBP_POSITIVE_REGULATION_OF_CAMP_PKA_SIGNAL_TRANSDUCTION',
        'GOBP_NEGATIVE_REGULATION_OF_CAMP_MEDIATED_SIGNALING',
        'GOBP_NEGATIVE_REGULATION_OF_CAMP_PKA_SIGNAL_TRANSDUCTION',
        'GOBP_CAMP_BIOSYNTHETIC_PROCESS',
        'GOBP_CAMP_CATABOLIC_PROCESS',
        'GOBP_CAMP_MEDIATED_SIGNALING',
        'GOBP_CAMP_METABOLIC_PROCESS',
        'GOBP_CAMP_PKA_SIGNAL_TRANSDUCTION',
        'GOBP_NEGATIVE_REGULATION_OF_CREB_TRANSCRIPTION_FACTOR_ACTIVITY',
        'GOBP_POSITIVE_REGULATION_OF_CREB_TRANSCRIPTION_FACTOR_ACTIVITY'
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
                for day in ['d10','d14','d21']:
                    if f'{day}_BI' in groups and f'{day}_Cont' in groups:
                        sc.tl.rank_genes_groups(da,groupby='sample_group',groups=[f'{day}_BI'],reference=f'{day}_Cont')
                        df = sc.get.rank_genes_groups_df(da,group=None)
                        df_sig = df[df.pvals_adj<0.05]
                        df_up = df_sig[(df_sig.logfoldchanges>0) & (df_sig.logfoldchanges<10)]
                        df_dw = df_sig[(df_sig.logfoldchanges<0) & (df_sig.logfoldchanges>-10)]
                        df_up.to_csv(join(oodir,rename(name)+'_'+day+'_up.csv'))
                        df_dw.to_csv(join(oodir,rename(name)+'_'+day+'_dw.csv'))
                        print(len(df_up))
                        print(len(df_dw))
                        if len(df_up) < 2 or  len(df_dw) < 2:
                            continue
                        #enr_up = gp.enrichr(df_up.names,gene_sets=gmt,outdir=None,organism='mouse',background=bg,verbose=True)
                        #enr_dw = gp.enrichr(df_dw.names,gene_sets=gmt,outdir=None,organism='mouse',background=bg,verbose=True)
                        enr_up = gp.enrichr(df_up.names,gene_sets={t:gmt[t] for t in target},outdir=None,organism='mouse',background=bg,verbose=True)
                        enr_dw = gp.enrichr(df_dw.names,gene_sets={t:gmt[t] for t in target},outdir=None,organism='mouse',background=bg,verbose=True)
                        print(enr_up.res2d)
                        print(enr_dw.res2d)
                        if enr_up.res2d is None or enr_dw.res2d is None:
                            continue
                        if len(enr_up.res2d) < 2 or len(enr_dw.res2d) < 2:
                            continue
                        enr_up.res2d['UP_DW'] = 'UP'
                        enr_dw.res2d['UP_DW'] = 'DOWN'
                        #enr_up.res2d = enr_up.res2d[enr_up.res2d['Combined Score'] < 100]
                        #enr_dw.res2d = enr_dw.res2d[enr_dw.res2d['Combined Score'] < 100]
                        enr_up.res2d = enr_up.res2d.sort_values('Combined Score',ascending=False)
                        enr_dw.res2d = enr_dw.res2d.sort_values('Combined Score',ascending=False)
                        enr_res = pd.concat([enr_up.res2d.head(10), enr_dw.res2d.head(10)])
                        #df.to_csv(join(oodir,rename(name)+'_'+day+'.csv'),index=False)
                        try:
                            enr_res.to_csv(join(oodir,rename(name)+'_'+day+'.csv'),index=False)
                            gp.dotplot(enr_up.res2d,figsize=(3,5),
                                cmap=NbDr,show_ring=True,top_term=20,size=10,cutoff=0.2,
                                title=f'{name} {day}_BI vs {day}_Cont',
                                ofname=join(oodir,rename(name)+'_'+day+'.png')
                            )
                        except ValueError:
                            pass

if __name__ == "__main__" :
    main()
