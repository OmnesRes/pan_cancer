##script for creating a histogram

## Load necessary modules
import pylab as plt
import numpy as np
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)

pvalues=map(float,pvalues)

##decide how man bins, 100 is the maximum possible due to only having two sig figs
number=100.0

counts={}

##use a dictionary to populate the bins
for i in range(int(number)):
    for j in pvalues:
        if i/number<j<=(i+1)/number:
            counts[i]=counts.get(i,0)+1

##convert the dictionary to a list
mylist=zip(counts.keys(),counts.values())

##sort the list so that the bins are in order
mylist.sort()



##plot the data with pylab
fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=.2)

ax.bar([i[0]/number for i in mylist],[i[1] for i in mylist],color='b',width=1/number,linewidth=2.0)
ax.set_xlim((0,1))

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
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xticks([i/10.0 for i in range(0,11)])
ax.set_xticklabels(['0']+[str(i/10.0) for i in range(1,11)])
ax.set_ylabel('Frequency',fontsize=60,labelpad=20)
ax.set_xlabel('Raw Cox P-value',fontsize=60,labelpad=20)
plt.show()
