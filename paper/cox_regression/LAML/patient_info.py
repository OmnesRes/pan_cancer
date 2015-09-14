## A script for extracting info about the patients used in the analysis

## Load necessary modules

from rpy2 import robjects as ro
import numpy as np
import os
ro.r('library(survival)')


##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

##LAML did not contain a follow_up file, which might explain why cbioportal doesn't allow kaplans for this cancer
##However there was a lot of clinical data in the clinical_patient file

## Sex and age information were taken from the "clinical_patient" file.  A dictionary was created for sex.
clinical=[]

more_clinical={}

sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1


f=open(os.path.join(BASE_DIR,'tcga_data','LAML','clinical','nationwidechildrens.org_clinical_patient_laml.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
for i in data:
    try:
        more_clinical[i[0]]=[0,sex_dict[i[4]],int(i[24])]
        if i[15]=='Alive':
            clinical.append([i[0],int(i[16]),'Alive'])
        elif i[15]=='Dead':
            clinical.append([i[0],int(i[17]),'Dead'])
        else:
            pass

    except:
        pass

## only patients who had a follow up time greater than 0 days are included in the analysis
clinical=[i for i in clinical if i[1]>0]


final_clinical=[]

## A new list containing both follow up times and sex and age is constructed.
## Only patients with sex and age information are included.
## Data is [[Patient ID, time (days), vital status, 0, sex, age at diagnosis],...]

for i in clinical:
    if i[0] in more_clinical:
        final_clinical.append(i+more_clinical[i[0]])



## Need to map the mRNA files to the correct patients
## The necessary information is included in the FILE_SAMPLE_MAP.txt file
f=open(os.path.join(BASE_DIR,'tcga_data','LAML','FILE_SAMPLE_MAP.txt'))
f.readline()
data=[i.strip().split() for i in f if i!='\n']



## 03 indicates a primary blood tumor, and only primary tumors are included in this analysis
TCGA_to_mrna={}
for i in data:
    ## The normalized data files are used
    if 'genes.normalized_results' in i[0]:
        if i[1].split('-')[3][:-1]=='03':
            x=''.join([k+j for k,j in zip(['','-','-'],i[1].split('-')[:3])])
            TCGA_to_mrna[x]=TCGA_to_mrna.get(x,[])+[i[0]]



clinical_and_files=[]
## We only care about patients that contained complete clinical information
for i in final_clinical:
    if TCGA_to_mrna.has_key(i[0]):
        ## The mRNA files are added to the clinical list
        ## Data structure: [[Patient ID, time (days), vital status, 0, sex, age at diagnosis,[mRNA files]],...]
        clinical_and_files.append(i+[TCGA_to_mrna[i[0]]])
    else:
        pass



##print average age at diagnosis
age=np.mean([i[5] for i in clinical_and_files])

##print number of males
males=len([i for i in clinical_and_files if i[4]==0])

##print number of females
females=len([i for i in clinical_and_files if i[4]==1])

##to get the median survival we need to call survfit from r


##prepare variables for R
ro.globalenv['times']=ro.IntVector([i[1] for i in clinical_and_files])

##need to create a dummy variable group
ro.globalenv['group']=ro.IntVector([0 for i in clinical_and_files])

##need a vector for deaths
death_dic={}
death_dic['Alive']=0
death_dic['Dead']=1
ro.globalenv['died']=ro.IntVector([death_dic[i[2]] for i in clinical_and_files])

res=ro.r('survfit(Surv(times,died) ~ as.factor(group))')

#the number of events(deaths) is the fourth column of the output
deaths=str(res).split('\n')[-2].strip().split()[3]


#the median survival time is the fifth column of the output
median=str(res).split('\n')[-2].strip().split()[4]

##write data to a file
f=open('patient_info.txt','w')
f.write('Average Age')
f.write('\t')
f.write('Males')
f.write('\t')
f.write('Females')
f.write('\t')
f.write('Deaths')
f.write('\t')
f.write('Median Survival')
f.write('\n')

f.write(str(age))
f.write('\t')
f.write(str(males))
f.write('\t')
f.write(str(females))
f.write('\t')
f.write(deaths)
f.write('\t')
f.write(median)







