##script for finding the overlap in the top 100 most significant gene sets from msigdb for good and bad genes

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




##get the 100 most enriched protective and harmful gene sets for each cancer
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','good_overlap'))
BLCA_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            BLCA_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','bad_overlap'))
BLCA_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            BLCA_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','LGG','good_overlap'))
LGG_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LGG_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','LGG','bad_overlap'))
LGG_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LGG_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','good_overlap'))
BRCA_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            BRCA_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','bad_overlap'))
BRCA_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            BRCA_bad.append(x.split('\t')[0])
            x=f.readline()

            

f=open(os.path.join(BASE_DIR,'cox_regression','CESC','good_overlap'))
CESC_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            CESC_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','CESC','bad_overlap'))
CESC_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            CESC_bad.append(x.split('\t')[0])
            x=f.readline()
            



f=open(os.path.join(BASE_DIR,'cox_regression','COAD','good_overlap'))
COAD_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            COAD_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','COAD','bad_overlap'))
COAD_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            COAD_bad.append(x.split('\t')[0])
            x=f.readline()
            

f=open(os.path.join(BASE_DIR,'cox_regression','GBM','good_overlap'))
GBM_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            GBM_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','GBM','bad_overlap'))
GBM_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            GBM_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','good_overlap'))
HNSC_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            HNSC_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','bad_overlap'))
HNSC_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            HNSC_bad.append(x.split('\t')[0])
            x=f.readline()
            

f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','good_overlap'))
KIRC_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            KIRC_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','bad_overlap'))
KIRC_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            KIRC_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','good_overlap'))
KIRP_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            KIRP_good.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','bad_overlap'))
KIRP_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            KIRP_bad.append(x.split('\t')[0])
            x=f.readline()


f=open(os.path.join(BASE_DIR,'cox_regression','LAML','good_overlap'))
LAML_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LAML_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','LAML','bad_overlap'))
LAML_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LAML_bad.append(x.split('\t')[0])
            x=f.readline()




f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','good_overlap'))
LIHC_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LIHC_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','bad_overlap'))
LIHC_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LIHC_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','good_overlap'))
LUAD_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LUAD_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','bad_overlap'))
LUAD_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LUAD_bad.append(x.split('\t')[0])
            x=f.readline()


f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','good_overlap'))
LUSC_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LUSC_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','bad_overlap'))
LUSC_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            LUSC_bad.append(x.split('\t')[0])
            x=f.readline()

f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','good_overlap'))
SKCM_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            SKCM_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','bad_overlap'))
SKCM_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            SKCM_bad.append(x.split('\t')[0])
            x=f.readline()


f=open(os.path.join(BASE_DIR,'cox_regression','OV','good_overlap'))
OV_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            OV_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','OV','bad_overlap'))
OV_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            OV_bad.append(x.split('\t')[0])
            x=f.readline()



f=open(os.path.join(BASE_DIR,'cox_regression','STAD','good_overlap'))
STAD_good=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            STAD_good.append(x.split('\t')[0])
            x=f.readline()
            
f=open(os.path.join(BASE_DIR,'cox_regression','STAD','bad_overlap'))
STAD_bad=[]
x=f.readline()
while x!='':
    x=f.readline()
    if x=='Gene Set Name\t# Genes in Gene Set (K)\tDescription\t# Genes in Overlap (k)\tk/K\tp-value\tFDR q-value\n':
        x=f.readline()
        while x!='\n':
            STAD_bad.append(x.split('\t')[0])
            x=f.readline()




all_cancers1=[BLCA_good,BRCA_good,CESC_good,COAD_good,GBM_good,\
             HNSC_good,KIRC_good,KIRP_good,LAML_good,LGG_good,LIHC_good,\
             LUAD_good,LUSC_good,OV_good,SKCM_good,STAD_good]


all_cancers2=[BLCA_bad,BRCA_bad,CESC_bad,COAD_bad,GBM_bad,\
             HNSC_bad,KIRC_bad,KIRP_bad,LAML_bad,LGG_bad,LIHC_bad,\
             LUAD_bad,LUSC_bad,OV_bad,SKCM_bad,STAD_bad]



##create a list of lists of the overlaps, use all_cancers1 for good overlaps, all_cancers2 for bad overlaps
final_array=[]
for i in all_cancers2[::-1]:
    temp=[]
    for j in all_cancers2[::-1]:
        temp.append(len([k for k in j if k in i]))
    final_array.append(temp)

    

##plotting, use blue_yellow_red1 cmap for good overlaps, blue_yellow_red2 for bad overlaps
blue_yellow_red1 = make_colormap({0:'#00005C',.05:'#0000D0',.14:'#01BBCF',.15:'#33CC33',.2:'#FFFF00',.27:'#FF9900',.33:'#B47603',.35:'#A32900',1:'#751E00'})
blue_yellow_red2 = make_colormap({0:'#00005C',.05:'#0000D0',.14:'#01BBCF',.15:'#33CC33',.25:'#FFFF00',.3:'#FF9900',.38:'#B47603',.45:'#A32900',1:'#751E00'})
Z=np.array(final_array)
mask=np.tri(Z.shape[0],k=-1)
Z= np.ma.array(Z, mask=mask)
fig = plt.figure()
fig.subplots_adjust(bottom=.15)
fig.subplots_adjust(left=.15)
ax = fig.add_subplot(111)
figure=ax.imshow(Z,cmap=blue_yellow_red2,interpolation="nearest")
cbar=fig.colorbar(figure,pad=.02)
cbar.ax.tick_params(labelsize=40)
cbar.set_label('number of genes', rotation=270,fontsize=80,labelpad=25)
ax.set_yticks([i for i in range(0,16)])
ax.set_yticklabels(['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD'][::-1])
ax.tick_params(axis='y',labelsize=40)
ax.set_xticks([i for i in range(0,16)])
ax.set_xticklabels(['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD'][::-1],rotation=90)
ax.tick_params(axis='x',labelsize=40,)
ax.tick_params(axis='x',length=0,width=0)
ax.tick_params(axis='y',length=0,width=0)
ax.invert_yaxis()
ax.invert_xaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.show()




##get the overlaps within cancers


final_array=[]
for index1,i in enumerate(all_cancers1):
    for index2,j in enumerate(all_cancers2):
        if index2==index1:
            final_array.append(len([k for k in j if k in i]))

f=open('overlap_within_cancers.txt','w')
f.write(str(final_array))
f.close()





