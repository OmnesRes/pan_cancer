##script for creating table_S2
##Need to find the most common good and bad gene sets from figure 2, and ignore non-specific gene sets

##load the necessary modules
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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



##we don't want any gene sets that are shared between the good and bad gene sets within a cancer
common_sets={}
for index1,i in enumerate(all_cancers1):
    for index2,j in enumerate(all_cancers2):
        if index2==index1:
            for k in j:
                if k in i:
                    common_sets[k]=''



names=['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD']



bad_counts={}
for index, i in enumerate(all_cancers2):
    for j in i:
        bad_counts[j]=bad_counts.get(j,[])+[names[index]]

good_counts={}
for index, i in enumerate(all_cancers1):
    for j in i:
        good_counts[j]=good_counts.get(j,[])+[names[index]]


badlist=[[len(i),j,k] for i,j,k in zip(bad_counts.values(),bad_counts.values(),bad_counts.keys())]
badlist.sort(reverse=True)


goodlist=[[len(i),j,k] for i,j,k in zip(good_counts.values(),good_counts.values(),good_counts.keys())]
goodlist.sort(reverse=True)

f=open('bad.txt','w')
for i in badlist:
    f.write(i[-1])
    f.write('\t')
    f.write(str(i[0]))
    f.write('\t')
    f.write(''.join(k+j for k,j in zip(i[1],[',']*len(i[1])))[:-1])
    f.write('\t')
    f.write(str(len(good_counts.get(i[-1],[]))))
    f.write('\t')
    f.write(''.join(k+j for k,j in zip(good_counts.get(i[-1],[]),[',']*len(good_counts.get(i[-1],[]))))[:-1])
    f.write('\n')

f.close()

f=open('good.txt','w')
for i in goodlist:
    f.write(i[-1])
    f.write('\t')
    f.write(str(i[0]))
    f.write('\t')
    f.write(''.join(k+j for k,j in zip(i[1],[',']*len(i[1])))[:-1])
    f.write('\t')
    f.write(str(len(bad_counts.get(i[-1],[]))))
    f.write('\t')
    f.write(''.join(k+j for k,j in zip(bad_counts.get(i[-1],[]),[',']*len(bad_counts.get(i[-1],[]))))[:-1])
    f.write('\n')

f.close()




f=open('common.txt','w')
for i in common_sets:
    f.write(i)
    f.write('\n')
f.close()
