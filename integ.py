import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    data = sc.read_h5ad('all.h5')
    data.obs_names_make_unique()
    sc.pp.scrublet(data, batch_key="sample_name")
    data.layers['counts'] = data.X.copy()
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=2000, batch_key="sample_name")
    sc.pl.highly_variable_genes(data,show=False,save='_hvg.png',)
    sc.tl.pca(data)
    sce.pp.harmony_integrate(data,'sample_name')
    
    sc.pp.neighbors(data,use_rep='X_pca_harmony')
    sc.tl.umap(data)
    sc.pl.umap(
        data,
        color="sample_name",
        size=2,
        show=False,
        save='_umap.png'
    )
    sc.tl.leiden(data, flavor="igraph", n_iterations=2)
    sc.pl.umap(
        data,
        color=["leiden", "predicted_doublet", "doublet_score"],
        wspace=0.5,
        size=3,
        show=False,
        save='_umap2.png'
    )
    sc.pl.umap(
        data,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
        show=False,
        save='_umap3.png'
    )
    ress = np.arange(0.1,2,0.1)
    for res in ress:
        sc.tl.leiden(
            data, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
        )
    sc.pl.umap(
        data,color=[f"leiden_res_{res:4.2f}" for res in ress],legend_loc="on_data",legend_fontsize='xx-small',ncols=5
    )
    data.write('all2.h5')

if __name__ == "__main__" :
    main()
