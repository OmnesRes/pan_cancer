import numpy as np
import math
f=open('coeffs_pvalues_adjusted.txt')
data=[i.strip().split() for i in f]
genes,coeffs,pvalues,adjusted=zip(*data)

for_u_pos=[float(i[1]) for i in data if float(i[1])>=0]
for_u_neg=[float(i[1]) for i in data if float(i[1])<=0]

u_pos=np.percentile(for_u_pos,95)
u_neg=np.percentile([-1*i for i in for_u_neg],95)


normalized=[]
for i in data:
    x=float(i[1])
    if x>0:
        value=(2/(1+math.e**(-2*x/u_pos)))-1
    elif x<0:
        value=(2/(1+math.e**(-2*x/u_neg)))-1
    else:
        value=0
    normalized.append(value)

f=open('coeffs_normalized_pvalues_adjusted.txt','w')
for i,j,k,l,m in zip(genes,coeffs,normalized,pvalues,adjusted):
    f.write(i)
    f.write('\t')
    f.write(j)
    f.write('\t')
    f.write(str(k))
    f.write('\t')
    f.write(l)
    f.write('\t')
    f.write(m)
    f.write('\n')
f.close()

