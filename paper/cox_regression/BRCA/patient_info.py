## A script for extracting info about the patients used in the analysis

## Load necessary modules

from rpy2 import robjects as ro
import numpy as np
import os
ro.r('library(survival)')


##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


## There were three clinical files with nonredundant data.  V4.0 was found to be most up to date.
## V2.1 was more up to date than V1.5.   All three files are loaded with the more up to date file getting preference.


f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
## A patient can be listed multiple times in the file. The most recent listing (furthest down in the file), contains the most recent
## follow up data.  This code checks if the patient has already been loaded into the list, and if so, takes the more recent data.
## This required an empty value in the list initialization.
## Data is: [[Patient ID, time(days), Vital status],[Patient ID, time(days), Vital status],...]
clinical3=[['','','']]
for i in data:
    if clinical3[-1][0]==i[0]:
        if i[8]=='Alive':
            clinical3[-1]=[i[0],int(i[9]),'Alive']
        elif i[8]=='Dead':
            try:
                clinical3[-1]=[i[0],int(i[10]),'Dead']
            except:
                pass
        else:
            pass
    else:
        if i[8]=='Alive':
            clinical3.append([i[0],int(i[9]),'Alive'])
        elif i[8]=='Dead':
            try:
                clinical3.append([i[0],int(i[10]),'Dead'])
            except:
                pass
        else:
            pass


## Removing the empty value.
clinical=clinical3[1:]





f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical1=[['','','']]
for i in data:
    if i[0] not in [j[0] for j in clinical]:
        if clinical1[-1][0]==i[0]:
            if i[6]=='Alive':
                clinical1[-1]=[i[0],int(i[7]),'Alive']
            elif i[6]=='Dead':
                try:
                    clinical1[-1]=[i[0],int(i[8]),'Dead']
                except:
                    pass
            else:
                pass
        else:
            if i[6]=='Alive':
                clinical1.append([i[0],int(i[7]),'Alive'])
            elif i[6]=='Dead':
                try:
                    clinical1.append([i[0],int(i[8]),'Dead'])
                except:
                    pass
            else:
                pass

##merging data and removing the empty value
clinical+=clinical1[1:]


f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical2=[['','','']]
for i in data:
    if i[0] not in [j[0] for j in clinical]:
        if clinical2[-1][0]==i[0]:
            try:
                if i[6]=='Alive':
                    clinical2[-1]=[i[0],int(i[7]),'Alive']
                elif i[6]=='Dead':
                    try:
                        clinical2[-1]=[i[0],int(i[8]),'Dead']
                    except:
                        pass
                else:
                    pass
            except:
                pass
        else:
            try:
                if i[6]=='Alive':
                    clinical2.append([i[0],int(i[7]),'Alive'])
                elif i[6]=='Dead':
                    try:
                        clinical2.append([i[0],int(i[8]),'Dead'])
                    except:
                        pass
                else:
                    pass
            except:
                pass


##merging data and removing the empty value
clinical+=clinical2[1:]



## Grade, sex, and age information were taken from the "clinical_patient" file.  A dictionary was created for sex and grade.
more_clinical={}
grade_dict={}
grade_dict['Infiltrating Ductal Carcinoma']=1
grade_dict['Metaplastic Carcinoma']=3
grade_dict['Mucinous Carcinoma']=4
grade_dict['Medullary Carcinoma']=5
grade_dict['Infiltrating Lobular Carcinoma']=6




sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1

## The "clinical_patient" file can also contain patients not listed in the follow_up files.
## In these cases the clinical data for these patients gets appended to a new clinical list.

f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_patient_brca.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical4=[]
for i in data:
    try:
        more_clinical[i[0]]=[grade_dict[i[-11]],sex_dict[i[6]],int(i[20])]
    except:
        pass
    if i[13]=='Alive':
        clinical4.append([i[0],int(i[14]),'Alive'])
    elif i[13]=='Dead':
        try:
            clinical4.append([i[0],int(i[15]),'Dead'])
        except:
            pass
    else:
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



## Need to map the mRNA files to the correct patients
## The necessary information is included in the FILE_SAMPLE_MAP.txt file
f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','FILE_SAMPLE_MAP.txt'))
f.readline()
data=[i.strip().split() for i in f if i!='\n']

## 01 indicates a primary tumor, and only primary tumors are included in this analysis
TCGA_to_mrna={}
for i in data:
    ## The normalized data files are used
    if 'genes.normalized_results' in i[0]:
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







