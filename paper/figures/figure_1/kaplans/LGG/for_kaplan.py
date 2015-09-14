## A script for obtaining the normalized expression values of genes of interest and preparing them for R

## Load necessary modules
import os


##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))


f=open(os.path.join(BASE_DIR,'tcga_data','LGG','clinical','nationwidechildrens.org_clinical_follow_up_v1.0_lgg.txt'))
f.readline()
f.readline()
f.readline()
data=[i.split('\t') for i in f]
## A patient can be listed multiple times in the file. The most recent listing (furthest down in the file), contains the most recent
## follow up data.  This code checks if the patient has already been loaded into the list, and if so, takes the more recent data.
## This required an empty value in the list initialization.
## Data is: [[Patient ID, time(days), Vital status],[Patient ID, time(days), Vital status],...]
clinical=[['','','']]
for i in data:
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

## Removing the empty value.
clinical=clinical[1:]

## Grade, sex, and age information were taken from the "clinical_patient" file.  A dictionary was created for sex and grade.
more_clinical={}
grade_dict={}
grade_dict['G2']=2
grade_dict['G3']=3


sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1

## The "clinical_patient" file can also contain patients not listed in the follow_up files.
## In these cases the clinical data for these patients gets appended to a new clinical list.

f=open(os.path.join(BASE_DIR,'tcga_data','LGG','clinical','nationwidechildrens.org_clinical_patient_lgg.txt'))
f.readline()
f.readline()
f.readline()
clinical4=[]
data=[i.split('\t') for i in f]
for i in data:
    more_clinical[i[0]]=[grade_dict[i[4]],sex_dict[i[10]],int(i[-12])]
    if i[39]=='Alive':
        clinical4.append([i[0],int(i[40]),'Alive'])
    elif i[39]=='Dead':
        clinical4.append([i[0],int(i[41]),'Dead'])
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


## A new list containing both follow up times and grade, sex, and age is constructed.
## Only patients with grade, sex, and age information are included.
## Data is [[Patient ID, time (days), vital status, grade, sex, age at diagnosis],...]


final_clinical=[]

for i in clinical:
    if i[0] in more_clinical:
        final_clinical.append(i+more_clinical[i[0]])




## Need to map the mRNA files to the correct patients
## The necessary information is included in the FILE_SAMPLE_MAP.txt file
f=open(os.path.join(BASE_DIR,'tcga_data','LGG','FILE_SAMPLE_MAP.txt'))
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


death_dic={}
death_dic['Alive']=0
death_dic['Dead']=1

##need to get the column order of the clustergram
f=open(os.path.join(BASE_DIR,'figures','figure_1','clustergrams','LGG','column_order.txt'))
##double loop list comprehension
new_order=[int(j) for i in f for j in i.strip().split()]

##get the original order of the patients
f=open(os.path.join(BASE_DIR,'figures','figure_1','clustergrams','LGG','for_clustering_100_good_100_bad.txt'))
original_order=f.readline().split()


##using the column order, rearrange the patients so that they are the order in the clustergram
order_dict={}
for index,i in enumerate(original_order):
    order_dict[index+1]=i

final_order=[order_dict[i] for i in new_order]

##split the groups at TCGA-DU-7019


finaldata=[]
for i in clinical_and_files:
    if i[0] in final_order[:final_order.index('TCGA-DU-7019')]:
        finaldata.append([i[1],i[2],'1'])
    elif i[0] in final_order[final_order.index('TCGA-DU-7019'):]:
        finaldata.append([i[1],i[2],'2'])
    else:
        pass


f=open('for_kaplan.txt','w')
f.write('time')
f.write('\t')
f.write('Died')
f.write('\t')
f.write('group')
f.write('\n')
for i in finaldata:
    f.write(str(i[0]))
    f.write('\t')
    f.write(str(death_dic[i[1]]))
    f.write('\t')
    f.write(i[2])
    f.write('\n')
f.close()


