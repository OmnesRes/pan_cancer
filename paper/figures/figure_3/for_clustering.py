##script for finding genes significant in multiple cancers and getting the normalized Cox coefficients

##load necessary modules
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

##will need the adjusted pvalues and the normalized Cox Coefficients for each cancer
f=open(os.path.join(BASE_DIR,'cox_regression','BLCA','coeffs_normalized_pvalues_adjusted.txt'))
BLCA=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','BRCA','coeffs_normalized_pvalues_adjusted.txt'))
BRCA=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','CESC','coeffs_normalized_pvalues_adjusted.txt'))
CESC=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','COAD','coeffs_normalized_pvalues_adjusted.txt'))
COAD=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','GBM','coeffs_normalized_pvalues_adjusted.txt'))
GBM=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','HNSC','coeffs_normalized_pvalues_adjusted.txt'))
HNSC=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','KIRC','coeffs_normalized_pvalues_adjusted.txt'))
KIRC=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','KIRP','coeffs_normalized_pvalues_adjusted.txt'))
KIRP=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','LAML','coeffs_normalized_pvalues_adjusted.txt'))
LAML=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','LGG','coeffs_normalized_pvalues_adjusted.txt'))
LGG=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','LIHC','coeffs_normalized_pvalues_adjusted.txt'))
LIHC=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','LUAD','coeffs_normalized_pvalues_adjusted.txt'))
LUAD=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','LUSC','coeffs_normalized_pvalues_adjusted.txt'))
LUSC=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','OV','coeffs_normalized_pvalues_adjusted.txt'))
OV=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','SKCM','coeffs_normalized_pvalues_adjusted.txt'))
SKCM=[i.strip().split() for i in f]

f=open(os.path.join(BASE_DIR,'cox_regression','STAD','coeffs_normalized_pvalues_adjusted.txt'))
STAD=[i.strip().split() for i in f]



all_cancers=[BLCA,BRCA,CESC,COAD,GBM,HNSC,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,OV,SKCM,STAD]

##use a dictionary to keep track of the significant genes
sig_genes={}
for index,i in enumerate(all_cancers):
    for j in i:
        if float(j[-1])<=.05:
            sig_genes[j[0]]=sig_genes.get(j[0],[])+[index]

            

##use a dictionary to keep track of the Cox coefficients
final_coeffs={}
for index,i in enumerate(all_cancers):
    for j in i:
        final_coeffs[j[0]]=final_coeffs.get(j[0],[])+[j[2]]



##write data to a file, avoid trailing tabs
f=open('for_clustering.txt','w')
f.write('\t')
for i in ['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','SKCM','STAD'][:-1]:
    f.write(i)
    f.write('\t')
f.write('STAD')
f.write('\n')


for i in final_coeffs:
    if i in sig_genes:
        if len(final_coeffs[i])==16 and len(sig_genes[i])>=4:
            f.write(i)
            f.write('\t')
            for j in final_coeffs[i][:-1]:
                f.write(j)
                f.write('\t')
            f.write(final_coeffs[i][-1])
            f.write('\n')

f.close()
        
             

