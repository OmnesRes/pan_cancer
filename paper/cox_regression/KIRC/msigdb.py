#this code sorts the genes by pvalues and prints out the 250 most significant good and bad genes
#these gene sets were then investigated at http://www.broadinstitute.org/gsea/msigdb/annotate.jsp
#C1,CGP,CP,CP:KEGG,MIR,TFT,CM,BP,C6 gene sets were checked
f=open('coeffs_pvalues.txt')
data=[i.strip().split() for i in f]

pvalues=[]
for i in data:
    pvalues.append([float(i[-1]),float(i[1]),i[0]])

pvalues.sort()

good=[i for i in pvalues if i[1]<0]
bad=[i for i in pvalues if i[1]>0]

for i in good[:250]:
    print i[2]

print
print
print '---------------------'
print
print


for i in bad[:250]:
    print i[2]
