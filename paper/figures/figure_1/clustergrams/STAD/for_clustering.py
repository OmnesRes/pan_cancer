## A script for obtaining the normalized expression values of genes of interest and preparing them for R

## Load necessary modules
import numpy as np
import os
from rpy2 import robjects as ro


##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

## Read the follow up data
## A patient can be listed multiple times in the file. The most recent listing (furthest down in the file), contains the most recent
## follow up data.  This code checks if the patient has already been loaded into the list, and if so, takes the more recent data.
## This required an empty value in the list initialization.
## Data is: [[Patient ID, time(days), Vital status],[Patient ID, time(days), Vital status],...]

f=open(os.path.join(BASE_DIR,'tcga_data','STAD','clinical','nationwidechildrens.org_clinical_follow_up_v1.0_stad.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical=[['','','']]
for i in data:
    try:
        if clinical[-1][0]==i[0]:
            if i[8]=='Alive':
                clinical[-1]=[i[0],int(i[9]),'Alive']
            elif i[8]=='Dead':
                clinical[-1]=[i[0],int(i[10]),'Dead']
            else:
                pass
        else:
            if i[8]=='Alive':
                clinical.append([i[0],int(i[9]),'Alive'])
            elif i[8]=='Dead':
                clinical.append([i[0],int(i[10]),'Dead'])
            else:
                pass
    except:
        pass


## Removing the empty value.
clinical=clinical[1:]


## Grade, sex and age information were taken from the "clinical_patient" file.  A dictionary was created for grade and sex.
more_clinical={}
grade_dict={}
grade_dict['G1']=1
grade_dict['G2']=2
grade_dict['G3']=3



sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1

## The "clinical_patient" file can also contain patients not listed in the follow_up files.
## In these cases the clinical data for these patients gets appended to a new clinical list.

f=open(os.path.join(BASE_DIR,'tcga_data','STAD','clinical','nationwidechildrens.org_clinical_patient_stad.txt'))
f.readline()
f.readline()
f.readline()
clinical4=[]
data=[i.split('\t') for i in f]
for i in data:
    try:
        more_clinical[i[0]]=[grade_dict[i[4]],sex_dict[i[7]],int(i[41])]
        if i[26]=='Alive':
            clinical4.append([i[0],int(i[27]),'Alive'])
        elif i[26]=='Dead':
            clinical4.append([i[0],int(i[28]),'Dead'])
        else:
            pass

    except:
        pass


new_clinical=[]


##It is possible that the clinical data in the clinical_patient file is more up to date than the follow_up files
##All the clinical data is merged checking which data is the most up to date
for i in clinical4:
    if i[0] not in [j[0] for j in clinical]:
        new_clinical.append(i)
    else:
        if i[1]<=clinical[[j[0] for j in clinical].index(i[0])][1]:
            new_clinical.append(clinical[[j[0] for j in clinical].index(i[0])])
        else:
            new_clinical.append(i)

##also do the reverse since clinical can contain patients not included in clinical4
for i in clinical:
    if i[0] not in [j[0] for j in new_clinical]:
        new_clinical.append(i)

## only patients who had a follow up time greater than 0 days are included in the analysis
clinical=[i for i in new_clinical if i[1]>0]


final_clinical=[]

## A new list containing both follow up times and grade, sex, and age is constructed.
## Only patients with grade, sex, and age information are included.
## Data is [[Patient ID, time (days), vital status, grade, sex, age at diagnosis],...]


for i in clinical:
    if i[0] in more_clinical:
        final_clinical.append(i+more_clinical[i[0]])


## 01 indicates a primary tumor, and only primary tumors are included in this analysis
f=open(os.path.join(BASE_DIR,'tcga_data','STAD','FILE_SAMPLE_MAP.txt'))
f.readline()
data=[i.strip().split() for i in f if i!='\n']


## 01 indicates a primary tumor, and only primary tumors are included in this analysis
TCGA_to_mrna={}
for i in data:
    ##the v2.gene.quantification files were used
    if 'v2.gene.quantification' in i[0]:
        if i[1].split('-')[3][:-1]=='01':
            x=''.join([k+j for k,j in zip(['','-','-'],i[1].split('-')[:3])])
            TCGA_to_mrna[x]=TCGA_to_mrna.get(x,[])+[i[0]]



clinical_and_files=[]
## We only care about patients that contained complete clinical information
for i in final_clinical:
    if TCGA_to_mrna.has_key(i[0]):
        ## The mRNA files are added to the clinical list
        ## Data structure: [[Patient ID, time (days), vital status, grade, sex, age at diagnosis,[mRNA files]],...]
        clinical_and_files.append(i+[TCGA_to_mrna[i[0]]])
    else:
        pass


##this assumes you ran cox_regression.py and saved final_genes
f=open(os.path.join(BASE_DIR,'cox_regression','STAD','final_genes.txt'))
final_genes=[eval(i.strip()) for i in f]


##final_genes is a list of gene expression values for each patient included in the study
##convert to a dictionary to be able to access values of interest
##the values are in the same order as the patients in clinical_and_files

final_dict={}
for i in final_genes:
    for j in i:
        final_dict[j[0]]=final_dict.get(j[0],[])+[j[1]]


f=open(os.path.join(BASE_DIR,'cox_regression','STAD','coeffs_pvalues.txt'))
data2=[i.strip().split() for i in f]

pvalues=[]
for i in data2:
    pvalues.append([float(i[-1]),float(i[1]),i[0]])


pvalues.sort()



##get the normalized expression values for the good genes
for_matrix=[]
x=0
for i in pvalues:
    if i[1]<0:
        ro.globalenv['expression']=ro.FloatVector(final_dict[i[-1]])
        res=ro.r('round(qnorm((rank(expression, na.last="keep")-0.5)/sum(!is.na(expression))), digit=5)')
        inverse_norm=list(res)
        for_matrix.append([i[-1]]+inverse_norm)
        x+=1
        if x==100:
            break


f=open('for_clustering_100_good_100_bad.txt','w')

##write the header (patient names)
f.write('\t')
##don't want a trailing tab, may cause a NA value in R
for i in clinical_and_files[:-1]:
    f.write(i[0])
    f.write('\t')
f.write(clinical_and_files[-1][0])
f.write('\n')


##write the rows for good genes
for i in for_matrix:
    ##don't want a trailing tab, may cause a NA value in R
    for j in i[:-1]:
        f.write(str(j))
        f.write('\t')
    f.write(str(i[-1]))
    f.write('\n')
f.close()


##need a file for a color bar
f=open('genes_with_prognosis.txt','w')
##just want the gene names
for i in for_matrix:
    f.write(i[0])
    f.write('\t')
f.close()

##get the normalized expression values for the bad genes
for_matrix=[]
x=0
for i in pvalues:
    if i[1]>0:
        ro.globalenv['expression']=ro.FloatVector(final_dict[i[-1]])
        res=ro.r('round(qnorm((rank(expression, na.last="keep")-0.5)/sum(!is.na(expression))), digit=5)')
        inverse_norm=list(res)
        for_matrix.append([i[-1]]+inverse_norm)
        x+=1
        if x==100:
            break

f=open('for_clustering_100_good_100_bad.txt','a')
##write the rows for bad genes
for i in for_matrix:
    ##don't want a trailing tab, may cause a NA value in R
    for j in i[:-1]:
        f.write(str(j))
        f.write('\t')
    f.write(str(i[-1]))
    f.write('\n')

f.close()


##finish writing color bar file
f=open('genes_with_prognosis.txt','a')
##just want the gene names
for i in for_matrix[:-1]:
    f.write(i[0])
    f.write('\t')
f.write(for_matrix[-1][0])
f.write('\n')
f.write(''.join([i+j for i,j in zip(['Good']*100,['\t']*100)]))
f.write(''.join([i+j for i,j in zip(['Bad']*99,['\t']*99)]))
f.write('Bad')
f.close()




