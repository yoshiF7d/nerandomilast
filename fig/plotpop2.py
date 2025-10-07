import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker

filein='all5.h5'
metafile='250828meta.tsv'
cluster='leiden_res_0.20'

order=[
    "Alveolar epithelial cell",
    "B cell",
    "Endothelial cell",
    "Fibroblast",
    "Macrophage/Dendric cell",
    "Mesothelial cell",
    "Monocyte",
    "Myofibroblast",
    "NK cell",
    "Neutrophil",
    "T cell",
]

data = ad.read_h5ad(filein)
meta = pd.read_csv(metafile,sep='\t',index_col=0)
data.obs[cluster] = meta[cluster]
data.obs['sample_group'] = meta['sample_group']
data = data[data.obs.sample_group!='d21_Nint']
data = data[~data.obs[cluster].str.startswith('Classical')]
data.obs.loc[data.obs[cluster].str.endswith('monocyte'),cluster] = 'Monocyte'
data = data[data.obs[cluster] != 'Unassigned']
data = data[data.obs.sample_group != 'Unassigned']
keys = np.sort(data.obs[cluster].unique())
colors = data.uns[f'{cluster}_colors']
cdic = dict(zip(keys,colors))

df = data.obs[[cluster,'sample_group']]
days = ['d0','d8','d10','d14']

tab = []
btab,ctab = [],[]

for day in days:
    d = df[df.sample_group.str.startswith(day)]
    p = d[cluster].value_counts()    
    if day == 'd8':
        pass
    if day in ['d10','d14']:
        db = d[d.sample_group.str.endswith('BI')]
        dc = d[d.sample_group.str.endswith('Cont')]
        pb = db[cluster].value_counts()
        pc = dc[cluster].value_counts()
        btab.append(pb[order].to_list())
        ctab.append(pc[order].to_list())

    print(p)
    tab.append(p[order].to_list())

tab = np.array(tab)
tab = tab.transpose()
btab = np.array(btab)
btab = btab.transpose()
ctab = np.array(ctab)
ctab = ctab.transpose()

#print(tab)
fig,ax = plt.subplots()
bot=np.zeros(2)
botb=np.zeros(2)
botc=np.zeros(2)
width = 0.35
for i,o in enumerate(reversed(order)):
    #print(tab[i,:2])
    #print(btab[i])
    p1 = ax.bar([0,1],tab[i,:2],width,label=o,color=cdic[o],bottom=bot)
    ax.bar(np.array([2,3])-width/2-0.05,btab[i],width,bottom=botb,color=p1[0].get_facecolor())
    ax.bar(np.array([2,3])+width/2+0.05,ctab[i],width,bottom=botc,color=p1[0].get_facecolor())
    botb += btab[i]
    botc += ctab[i]
    bot += tab[i,:2]

#formatter = mticker.EngFormatter()
#ax.yaxis.set_major_formatter(formatter)
ax.set_ylabel('Cell counts')
#ax.set_title('Normalized cell counts')
ax.set_xticks([0,1,2-width/2-0.05,2+width/2+0.05,3-width/2-0.05,3+width/2+0.05])
ax.set_xticklabels(['Day 0','Day 8','Day 10 BI','Day 10 Cont','Day 14 BI','Day 14 Cont'],rotation=-45,ha='left',rotation_mode='anchor')
#ax.legend(loc='upper left',bbox_to_anchor=(1,1))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1,1), loc='upper left')

fig.savefig('allcellcounts.png',bbox_inches='tight',dpi=300)

fig,ax = plt.subplots()
bot=np.zeros(4)
width = 0.35
for i,o in enumerate(order):
    ax.bar([0,1,2,3],tab[i],width,label=o,bottom=bot)
    bot += tab[i]

ax.set_ylabel('Normalized cell counts')
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['Day 0','Day 8','Day 10','Day 14'])
ax.legend(loc='upper left',bbox_to_anchor=(1,1))
fig.savefig('allcellcountsall.png',bbox_inches='tight')

