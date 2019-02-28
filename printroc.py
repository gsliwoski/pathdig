import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys

try:
    infile = sys.argv[1]
    predcol1 = sys.argv[2]
    predcol2 = sys.argv[3]
except:
    sys.exit("quickroc.py scorefile.csv predcol1 predcol2")

df = pd.read_csv(infile)
df = df[(~df['label'].isnull())]
with open("badstructs.ls") as infile:
    badstructs = [x.strip() for x in infile.readlines()]
df = df[~df.structid.isin(badstructs)]
df1 = df[~df[predcol1].isnull()]
df2 = df[~df[predcol2].isnull()]

#df1 = df1[df1[predcol1]!=0]
#df2 = df2[df2[predcol2]!=0]

#df = df[df.pathprox>0.0]
x = list()
y = list()
pred1 = list(df1.sort_values(predcol1,ascending=False)['label'])
pred2 = list(df2.sort_values(predcol2,ascending=False)['label'])
vid = list(df1.sort_values(predcol1,ascending=False)['variantID'])
sid = list(df1.sort_values(predcol1,ascending=False)['structid'])
pos1 = df1[df1['label']=='clinvar'].shape[0]
neg1 = df1.shape[0] - pos1
pos2 = df2[df2['label']=='clinvar'].shape[0]
neg2 = df2.shape[0] - pos2
tp1 = 0
fp1 = 0
tp2 = 0
fp2 = 0
n1 = df1.shape[0]
n2 = df2.shape[0]

plt.gca().set_aspect('equal',adjustable='box')
x1 = list()
y1 = list()
x2= list()
y2 = list()

for i in range(max(n1,n2)):
    try:
        if pred1[i] == 'clinvar':
            tp1 += 1
        else:
            fp1 += 1
        x1.append(float(fp1)/neg1)
        y1.append(float(tp1)/pos1)
        print "{} {} {} {} {}".format(float(fp1)/neg1,float(tp1)/pos1,pred1[i],vid[i],sid[i])
    except IndexError:
        pass
    try:        
        if pred2[i] == 'clinvar':
            tp2 += 1
        else:
            fp2 += 1
        x2.append(float(fp2)/neg2)
        y2.append(float(tp2)/pos2)
    except IndexError:
        pass            
   
plt.plot(x1,y1,color="blue")
plt.plot(x2,y2,color="orange")
patch1 = mpatches.Patch(color="blue",label=predcol1)
patch2 = mpatches.Patch(color="orange",label=predcol2)
plt.legend(handles=[patch1,patch2])
mn = max(n1,n2)
plt.plot([x/float(mn) for x in range(mn)],[x/float(mn) for x in range(mn)],color="red")
#plt.savefig("{}_{}_rocq2.png".format(predcol1,predcol2))
