##A script for creating a table


## Load necessary modules
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


##load data for each cancer, find total genes analyzed and significant genes, get patient info
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BLCA=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BRCA=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','CESC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','CESC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
CESC=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','COAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','COAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
COAD=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','GBM','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','GBM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
GBM=[genes_analyzed,sig]+data



f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
HNSC=[genes_analyzed,sig]+data



f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRC=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRP=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','LAML','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','LAML','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LAML=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','LGG','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','LGG','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LGG=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LIHC=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUAD=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUSC=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','OV','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','OV','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
OV=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
SKCM=[genes_analyzed,sig]+data


f=open(os.path.join(BASE_DIR,'cox_regression','STAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
genes_analyzed=len(ids)
sig=len([i for i in sorted(map(float,adjusted)) if i <=.05])
f=open(os.path.join(BASE_DIR,'cox_regression','STAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
STAD=[genes_analyzed,sig]+data




all_cancers=[BLCA,BRCA,CESC,COAD,GBM,HNSC,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,SKCM,OV,STAD]

names=['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','SKCM','OV','STAD']

f=open('table_1.txt','w')
for i,j in zip(all_cancers,names):
    f.write(j)
    f.write('\t')
    ##write total patients (add males and females)
    f.write(str(int(i[3])+int(i[4])))
    f.write('\t')
    ##write median survival
    f.write(i[6])
    f.write('\t')
    ##write events
    f.write(i[5])
    f.write('\t')
    ##write average age at diagnosis
    f.write(i[2])
    f.write('\t')
    ##write male/female
    f.write(i[3]+'/'+i[4])
    f.write('\t')
    ##write genes analyzed
    f.write(str(i[0]))
    f.write('\t')
    ##write number of sig genes
    f.write(str(i[1]))
    f.write('\n')

f.close()







