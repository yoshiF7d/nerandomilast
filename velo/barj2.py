import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

df = pd.read_csv('plotj2.out',sep='\t',header=None,names=['Sample','Value'])
df['Group'] = df['Sample'].str.rsplit('_',n=1).str[0]
df[['Day','Condition']] = df['Group'].str.split('_',expand=True)

summary_df = df.groupby(['Day','Condition'])['Value'].agg(['mean','std']).reset_index()
plot_df = summary_df.pivot(index='Day',columns='Condition',values=['mean','std'])
summary_df.to_csv('barplot2.csv',float_format='%.2e',index=None)

labels = plot_df.index
bi_means = plot_df[('mean','BI')]
cont_means = plot_df[('mean','Cont')]
bi_std = plot_df[('std','BI')]
cont_std = plot_df[('std','Cont')]

x = np.arange(len(labels))
width = 0.45

fig,ax = plt.subplots(figsize=(5,4))

#rects1 = ax.bar(x-width/2,bi_means,width,yerr=bi_std,label='BI',capsize=5,color='skyblue')
#rects2 = ax.bar(x+width/2,cont_means,width,yerr=cont_std,label='Cont',capsize=6,color='salmon')
rects1 = ax.bar(x-width/2,cont_means,width-0.1,yerr=cont_std/2,capsize=5,edgecolor='black',facecolor='#f94040')
rects2 = ax.bar(x+width/2,bi_means,width-0.1,yerr=bi_std/2,capsize=5,edgecolor='black',facecolor='#0000ff')
ax.scatter(x[0]-width/2,cont_means[0],marker="o",label='Day 10',color='black')
ax.scatter(x[0]+width/2,bi_means[0],marker="o",color='black')
ax.scatter(x[1]-width/2,cont_means[1],marker="s",label='Day 14',color='black')
ax.scatter(x[1]+width/2,bi_means[1],marker="s",color='black')
ax.scatter(x[2]-width/2,cont_means[2],marker="^",label='Day 21',color='black')
ax.scatter(x[2]+width/2,bi_means[2],marker="^",color='black')

ax.set_xticks(np.concatenate([x-width/2,x+width/2]))
ax.set_xticklabels(['placebo','placebo','placebo','Nerandomilast','Nerandomilast','Nerandomilast'],rotation=60,ha='right',va='center',rotation_mode='anchor')
ax.set_ylim([0,0.06])
ax.legend()
ax.axhline(0,color='gray',linewidth=0.8)

fig.tight_layout()
fig.savefig('barplot2.png',dpi=300)
