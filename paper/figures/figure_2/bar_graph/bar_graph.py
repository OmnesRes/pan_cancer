##script for creating a bar graph

##load necessary modules
import pylab as plt
import os


BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

f=open(os.path.join(BASE_DIR,'figures','figure_2','gene_set_overlap','overlap_within_cancers.txt'))

data=eval(f.read())

names=['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD']



fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=.2)

ax.bar(range(1,17),data,color='b',width=.9,linewidth=2.0,align='center')
ax.set_ylim(0,100)
ax.set_xlim((0.5,16.5))

for item in ax.get_yticklabels():
        item.set_fontsize(30)

for item in ax.get_xticklabels():
        item.set_fontsize(30)

ax.tick_params(axis='x',length=15,width=3,direction='out',labelsize=30)
ax.tick_params(axis='y',length=15,width=3,direction='out',labelsize=30)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['bottom'].set_position(['outward',10])
ax.spines['left'].set_position(['outward',10])
ax.spines['bottom'].set_bounds(1,16)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xticks(range(1,17))
ax.set_xticklabels(names,rotation=90)
ax.set_ylabel('Number of gene sets',fontsize=40,labelpad=20)
plt.show()
