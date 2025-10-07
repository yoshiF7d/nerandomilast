import anndata as ad
import pandas as pd
import scanpy as sc

filein='epi.h5ad'
cluster="Epithelium sub_resolution0.2"
data = ad.read_h5ad(filein)

data.obs['mycluster'] = data.obs[cluster].astype(str)
data.obs['mycluster'] = data.obs['mycluster'].map({'AT1':'AT1 cell','AT2':'AT2 cell','CLDN4+ cell':'PATS-1 cell','CLDN18+ cell':'PATS-2 cell'})
cluster = 'mycluster'

markers = {
    #'AT1':['Cav1','Ager','Aqp5','Hopx'],
    #'Transitional\n(CLDN4+, CLDN18+cell)':['Krt8','Krt19','Cldn4','Cldn18','Lgals3','Nupr1','Ddit3'],
    'AT1 cell':['Cav1','Ager','Aqp5','Pdpn'],
    'Transitional cell\n(PATS-1 cell, PATS-2 cell)':['Krt8','Krt19','Cldn4','Gclc','Malat1','Sfn','Cldn18','Prdx6','Clic3','Clic5'],
    'AT2 cell':['Sftpc','Sftpa1','Lamp3','Abca3'],
}

sc.pl.dotplot(data,var_names=markers,groupby=cluster,show=False,save='_dot.png')
sc.pl.heatmap(data,var_names=markers,groupby=cluster,show=False,save='_heat.png')
