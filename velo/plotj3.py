import argparse,tomllib
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import patches

colors = {
    'BI': {
        'd10': 'skyblue',
        'd14': 'dodgerblue',
        'd21': 'mediumblue'
    },
    'Cont': {
        'd10': 'sandybrown',
        'd14': 'darkorange',
        'd21': 'orangered'
    }
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--params',required=True)
    args = parser.parse_args()
    
    with open(args.params,'rb') as f:
        params=tomllib.load(f)

    rad = np.deg2rad(params['t'])
    R = np.array([[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]])
    
    fig,ax = plt.subplots()
    for day in ['d10','d14','d21']:
        for tre in ['BI','Cont']:
            name = f'{day}_{tre}'
            ps = []
            vs = []
            for rep in ['1','2']:
                data = ad.read_h5ad(f'result/{name}_{rep}/{name}_{rep}.h5ad')
                ps.append(data.obsm['X_umap'])
                vs.append(data.obsm['velocity_umap'])
            
            ps = np.concatenate(ps)
            vs = np.concatenate(vs)
            x,y = (R @ [ps[:,0]-params['x'],ps[:,1]-params['y']])
            ind = np.where((y/(0.5*params['a']))**2 + (x/(0.5*params['b']))**2 <= 1.0)
            jj = vs[ind,0]*np.cos(rad) - vs[ind,1]*np.sin(rad)
            jj = jj[0]
            
            ax.hist(jj,32,range=[-0.1,0.1],label=name,alpha=0.5,color=colors[tre][day])
            j = np.mean(jj)
            print('\t'.join([name,str(j),str(np.std(jj))]))
    
    ax.legend()
    fig.savefig(f'hist.png',dpi=200)

if __name__ == '__main__':
    main()
