##script for finding the overlap in the top 100 most significant genes in each cancer and plotting results

##load necessary modules
import pylab as plt
import numpy as np
import math
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


##I did not write this function, from http://depts.washington.edu/clawpack/clawpack-4.6.3/python/pyclaw/plotters/colormaps.py
##-------------------------
def make_colormap(colors):
##-------------------------
    """
    Define a new color map based on values specified in the dictionary
    colors, where colors[z] is the color that value z should be mapped to,
    with linear interpolation between the given values of z.

    The z values (dictionary keys) are real numbers and the values
    colors[z] can be either an RGB list, e.g. [1,0,0] for red, or an
    html hex string, e.g. "#ff0000" for red.
    """

    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort
    
    z = sort(colors.keys())
    n = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)
    
    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        #i'th color at level z[i]:
        Ci = colors[z[i]]      
        if type(Ci) == str:
            # a hex string of form '#ff0000' for example (for red)
            RGB = CC.to_rgb(Ci)
        else:
            # assume it's an RGB triple already:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])

    cmap_dict = {}
    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap



def compare3(first,second):
    if float(first[-1])>float(second[-1]):
        return 1
    elif float(first[-1])<float(second[-1]):
        return -1
    else:
        return 0


##get the 100 most significant genes for each cancer
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','coeffs_pvalues.txt'))
BLCA=[i.strip().split() for i in f]
BLCA.sort(cmp=compare3)
BLCA_dict_100={}
for i in BLCA[:100]:
    BLCA_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','LGG','coeffs_pvalues.txt'))
LGG=[i.strip().split() for i in f]
LGG.sort(cmp=compare3)
LGG_dict_100={}
for i in LGG[:100]:
    LGG_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','coeffs_pvalues.txt'))
BRCA=[i.strip().split() for i in f]
BRCA.sort(cmp=compare3)
BRCA_dict_100={}
for i in BRCA[:100]:
    BRCA_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','CESC','coeffs_pvalues.txt'))
CESC=[i.strip().split() for i in f]
CESC.sort(cmp=compare3)
CESC_dict_100={}
for i in CESC[:100]:
    CESC_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','COAD','coeffs_pvalues.txt'))
COAD=[i.strip().split() for i in f]
COAD.sort(cmp=compare3)
COAD_dict_100={}
for i in COAD[:100]:
    COAD_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','GBM','coeffs_pvalues.txt'))
GBM=[i.strip().split() for i in f]
GBM.sort(cmp=compare3)
GBM_dict_100={}
for i in GBM[:100]:
    GBM_dict_100[i[0]]=''

f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','coeffs_pvalues.txt'))
HNSC=[i.strip().split() for i in f]
HNSC.sort(cmp=compare3)
HNSC_dict_100={}
for i in HNSC[:100]:
    HNSC_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','coeffs_pvalues.txt'))
KIRC=[i.strip().split() for i in f]
KIRC.sort(cmp=compare3)
KIRC_dict_100={}
for i in KIRC[:100]:
    KIRC_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','coeffs_pvalues.txt'))
KIRP=[i.strip().split() for i in f]
KIRP.sort(cmp=compare3)
KIRP_dict_100={}
for i in KIRP[:100]:
    KIRP_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','LAML','coeffs_pvalues.txt'))
LAML=[i.strip().split() for i in f]
LAML.sort(cmp=compare3)
LAML_dict_100={}
for i in LAML[:100]:
    LAML_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','coeffs_pvalues.txt'))
LIHC=[i.strip().split() for i in f]
LIHC.sort(cmp=compare3)
LIHC_dict_100={}
for i in LIHC[:100]:
    LIHC_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','coeffs_pvalues.txt'))
LUAD=[i.strip().split() for i in f]
LUAD.sort(cmp=compare3)
LUAD_dict_100={}
for i in LUAD[:100]:
    LUAD_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','coeffs_pvalues.txt'))
LUSC=[i.strip().split() for i in f]
LUSC.sort(cmp=compare3)
LUSC_dict_100={}
for i in LUSC[:100]:
    LUSC_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','coeffs_pvalues.txt'))
SKCM=[i.strip().split() for i in f]
SKCM.sort(cmp=compare3)
SKCM_dict_100={}
for i in SKCM[:100]:
    SKCM_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','OV','coeffs_pvalues.txt'))
OV=[i.strip().split() for i in f]
OV.sort(cmp=compare3)
OV_dict_100={}
for i in OV[:100]:
    OV_dict_100[i[0]]=''
    
f=open(os.path.join(BASE_DIR,'cox_regression','STAD','coeffs_pvalues.txt'))
STAD=[i.strip().split() for i in f]
STAD.sort(cmp=compare3)
STAD_dict_100={}
for i in STAD[:100]:
    STAD_dict_100[i[0]]=''


all_cancers=[BLCA_dict_100,BRCA_dict_100,CESC_dict_100,COAD_dict_100,\
             GBM_dict_100,HNSC_dict_100,KIRC_dict_100,KIRP_dict_100,LAML_dict_100,\
             LGG_dict_100,LIHC_dict_100,LUAD_dict_100,LUSC_dict_100,OV_dict_100,\
             SKCM_dict_100,STAD_dict_100]

final_array=[]
for i in all_cancers[::-1]:
    temp=[]
    for j in all_cancers[::-1]:
        ##compute overlap
        temp.append(len([k for k in j if k in i]))
    final_array.append(temp)

##create a custom colormap
blue_yellow_red = make_colormap({0:'w',.05:'#85A3E0',.1:'#3366CC',.2:'#00FF00',.3:'#FFFF66',0.4:'#FF9966', 1:'#CC3300'})

##plot
Z=np.array(final_array)
mask=np.tri(Z.shape[0],k=-1)
Z= np.ma.array(Z, mask=mask)
fig = plt.figure()
fig.subplots_adjust(bottom=.15)
fig.subplots_adjust(left=.15)
ax = fig.add_subplot(111)
figure=ax.imshow(Z,cmap=blue_yellow_red,interpolation="nearest")
cbar=fig.colorbar(figure,pad=.02)
cbar.ax.tick_params(labelsize=40)
cbar.set_label('number of genes', rotation=270,fontsize=80,labelpad=25)
ax.set_yticks([i for i in range(0,16)])
ax.set_yticklabels(['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD'][::-1])
ax.tick_params(axis='y',labelsize=40)
ax.set_xticks([i for i in range(0,16)])
ax.set_xticklabels(['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD'][::-1],rotation=90)
ax.tick_params(axis='x',labelsize=40)
ax.tick_params(axis='x',length=0,width=0)
ax.tick_params(axis='y',length=0,width=0)
ax.invert_yaxis()
ax.invert_xaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.show()









