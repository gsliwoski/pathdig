import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import matplotlib.patches as mpatches
from sklearn import metrics

df = pd.read_csv("tmp3.csv",sep=" ")
cols = ['red','orange','green','blue','violet']
df['cat'] = pd.cut(df.nvar,5)
labels = list()
for i,x in enumerate(df.cat.unique()):
    plt.scatter(df[df.cat==x].fpr,df[df.cat==x].tpr,color=cols[i],s=70/(2*i+1))
    labels.append(mpatches.Patch(color=cols[i],label=str(x)))
plt.legend(handles=labels)

df2 = pd.read_csv("score_compare.csv")
df2 = df2[~df2.label.isnull()]
with open("badstructs.ls") as infile:
    badstructs = [x.strip() for x in infile.readlines()]
df2 = df2[~df2.structid.isin(badstructs)]
df2 = df2[~df2['pp_me'].isnull()]
y = df2['pp_me']
l = df2.label
fpr,tpr,t = metrics.roc_curve(l,y,pos_label="clinvar")
plt.plot(fpr,tpr,color="cyan",linewidth=.5)

plt.savefig("colorpp.png")


