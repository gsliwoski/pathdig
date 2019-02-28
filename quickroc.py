from sklearn import metrics
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    infile = sys.argv[1]
    predcol1 = sys.argv[2]
except:
    sys.exit("quickroc3.py infile predcol1 [predcol2]")

try:
    predcol2 = sys.argv[3]
except:
    predcol2 = predcol1

df = pd.read_csv(infile)
df = df[~df.label.isnull()]
try:
    with open("badstructs.ls") as infile:
        badstructs = [x.strip() for x in infile.readlines()]
except OSError:
    badstructs = list()
df = df[~df.structid.isin(badstructs)]
df1 = df[~df[predcol1].isnull()]
df2 = df[~df[predcol2].isnull()]

y1 = df1[predcol1]
y2 = df2[predcol2]
l1 = df1.label
l2 = df2.label

#print y1
#print l1
fpr1,tpr1,t1 = metrics.roc_curve(l1,y1,pos_label="clinvar")
fpr2,tpr2,t2 = metrics.roc_curve(l2,y2,pos_label="clinvar")
labels = list()            
plt.gca().set_aspect('equal',adjustable='box')
plt.plot(fpr1,tpr1,color="blue")
labels.append(mpatches.Patch(color="blue",label=predcol1))
if predcol1!=predcol2:
    plt.plot(fpr2,tpr2,color="orange")
    labels.append(mpatches.Patch(color="orange",label=predcol2))
mn = max(df1.shape[0],df2.shape[0])
plt.plot([x/float(mn) for x in range(mn)],[x/float(mn) for x in range(mn)],color="red")
plt.legend(handles=labels)
of = "{}_{}_roc.png".format(predcol1,predcol2) if predcol1!=predcol2 else "{}_roc.png".format(predcol1)
plt.savefig(of)

