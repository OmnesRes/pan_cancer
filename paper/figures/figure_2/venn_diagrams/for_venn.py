##script for parsing the msigdb outputs and writing results to a file for R

##load necessary modules

import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


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
            
f=open('COAD_good.txt','w')
for i in COAD_good:
    f.write(i)
    f.write('\n')
f.close()
            

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

f=open('GBM_good.txt','w')
for i in GBM_good:
    f.write(i)
    f.write('\n')
f.close()


            

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

f=open('KIRP_bad.txt','w')
for i in KIRP_bad:
    f.write(i)
    f.write('\n')
f.close()

            
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

f=open('LIHC_bad.txt','w')
for i in LIHC_bad:
    f.write(i)
    f.write('\n')
f.close()


            
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

f=open('LUAD_bad.txt','w')
for i in LUAD_bad:
    f.write(i)
    f.write('\n')
f.close()



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

f=open('LUSC_good.txt','w')
for i in LUSC_good:
    f.write(i)
    f.write('\n')
f.close()


