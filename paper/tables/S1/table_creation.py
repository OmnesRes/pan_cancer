##A script for creating tables for each cancer, with the data sorted 

def compare(first,second):
    if float(first[-2])>float(second[-2]):
        return 1
    elif float(first[-2])<float(second[-2]):
        return -1
    else:
        return 0

## Load necessary modules
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

##need to get the gene ids from a RNA-SEQV2 file, any file will work
f=open(os.path.join(BASE_DIR,'tcga_data','GBM','mrna','unc.edu.0cbec58e-f95e-4c60-a85d-210dc56bdf3c.1545137.rsem.genes.normalized_results'))
f.readline()
id_to_gene={}


data=[i.split()[0] for i in f]
for i in data:
    id_to_gene[i.split('|')[1]]=i.split('|')[0]


##load the data that will be in the table for each cancer and add ids
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
BLCA=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
BRCA=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','CESC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
CESC=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','COAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
COAD=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','GBM','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
GBM=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
HNSC=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
KIRC=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
KIRP=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','LAML','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
LAML=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','LGG','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
LGG=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
LIHC=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
LUAD=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
LUSC=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','OV','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
OV=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
SKCM=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)

f=open(os.path.join(BASE_DIR,'cox_regression','STAD','coeffs_normalized_pvalues_adjusted.txt'))
data=[i.strip().split() for i in f]
ids,coeffs,normalized,pvalues,adjusted=zip(*data)
STAD=zip(ids,[id_to_gene[i] for i in ids],coeffs,normalized,pvalues,adjusted)




all_cancers=[BLCA,LGG,BRCA,CESC,COAD,GBM,HNSC,KIRC,KIRP,LAML,LIHC,LUAD,LUSC,SKCM,OV,STAD]

names=['BLCA','LGG','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LIHC','LUAD','LUSC','SKCM','OV','STAD']

for i,j in zip(names,all_cancers):
    f=open(i+'.txt','w')
    for k in sorted(j,cmp=compare):
        for l in k:
            f.write(l)
            f.write('\t')
        f.write('\n')
    f.close()









