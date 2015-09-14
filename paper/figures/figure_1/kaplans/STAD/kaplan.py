import pylab as plt
from rpy2 import robjects as ro
ro.r('library(survival)')


##explanations of this code can be found at https://www.youtube.com/watch?v=cPEXL1YuhAk
f=open('for_kaplan.txt')
f.readline()
data=[i.strip().split() for i in f]
survtimes=[[int(i[0]),int(i[1])] for i in data if i[-1]=='1']

survtimes.sort(key=lambda x:(x[0],x[1]))


def kaplan(survtimes):
    h_coords=[]
    v_coords=[]
    lost=[]
    y=1
    x=0
    for i in survtimes:
        if i[1]!=1:
            lost.append([i[0],y])
        else:
            h_coords.append([i[0],y])
            y=1*len(survtimes[survtimes.index(i)+1:])/float(len(survtimes[survtimes.index(i):]))
            v_coords.append([i[0],h_coords[-1][-1],y])
            break
    newsurv=survtimes[survtimes.index(i)+1:]
    while len(newsurv)>0:
        newsurv,y,h_coords,v_coords,lost=loop(newsurv,y,h_coords,v_coords,lost)
    return (h_coords,v_coords,lost)



def loop(newsurv,y,h_coords,v_coords,lost):
    for j in newsurv:
        if j[1]!=1:
            lost.append([j[0],y])
        else:
            h_coords.append([j[0],y])
            y=y*len(newsurv[newsurv.index(j)+1:])/float(len(newsurv[newsurv.index(j):]))
            v_coords.append([j[0],h_coords[-1][-1],y])
            break
    newsurv=newsurv[newsurv.index(j)+1:]
    return (newsurv,y,h_coords,v_coords,lost)


k_plot=kaplan(survtimes)



fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=.1)
fig.subplots_adjust(top=.95)
fig.subplots_adjust(left=.05)
fig.subplots_adjust(right=.7)


width=3
start=0
for i in k_plot[0]:
    ax.hlines(i[1]*100,start,i[0],linewidths=width,color='b',label='cluster 1')
    start=i[0]

if k_plot[-1][-1][0]>k_plot[0][-1][0]:
    ax.hlines(k_plot[-1][-1][1]*100,k_plot[0][-1][0],k_plot[-1][-1][0],linewidths=width,color='b')

for i in k_plot[1]:
    cluster1=ax.vlines(i[0],i[2]*100-(width*67/1000.0),i[1]*100+(width*67/1000.0),linewidths=width,color='b')


for i in k_plot[2]:
    ax.vlines(i[0],(i[1]-.01)*100,(i[1]+.01)*100,color='b')




f=open('for_kaplan.txt')
f.readline()
data=[i.strip().split() for i in f]
survtimes=[[int(i[0]),int(i[1])] for i in data if i[-1]=='2']

survtimes.sort(key=lambda x:(x[0],x[1]))

k_plot=kaplan(survtimes)
start=0
for i in k_plot[0]:
    ax.hlines(i[1]*100,start,i[0],linewidths=width,color='r',label='cluster 2')
    start=i[0]

if k_plot[-1][-1][0]>k_plot[0][-1][0]:
    ax.hlines(k_plot[-1][-1][1]*100,k_plot[0][-1][0],k_plot[-1][-1][0],linewidths=width,color='r')

for i in k_plot[1]:
    cluster2=ax.vlines(i[0],i[2]*100-(width*67/1000.0),i[1]*100+(width*67/1000.0),linewidths=width,color='r')


for i in k_plot[2]:
    ax.vlines(i[0],(i[1]-.01)*100,(i[1]+.01)*100,color='r')





ax.tick_params(axis='x',length=15,width=3,direction='out',labelsize=30)
ax.tick_params(axis='y',length=15,width=3,direction='out',labelsize=30)
ax.spines['bottom'].set_position(['outward',10])
ax.spines['left'].set_position(['outward',10])
ax.set_xticks([0,800,1600,2400,3200,4000])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_bounds(0,100)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

ax.legend((cluster1,cluster2),('cluster 1','cluster 2'),frameon=False)

plt.ylim(0,105)
plt.xlim(0,)


plt.show()

##get the logrank p-value

f=open('for_kaplan.txt')
f.readline()
data=[i.strip().split() for i in f]
times=[int(i[0]) for i in data]
died=[int(i[1]) for i in data]
group=[int(i[2]) for i in data]

ro.globalenv['times']=ro.IntVector(times)
ro.globalenv['died']=ro.IntVector(died)
ro.globalenv['group']=ro.IntVector(group)

res=ro.r('survdiff(formula = Surv(times, died) ~ as.factor(group))')
print res











