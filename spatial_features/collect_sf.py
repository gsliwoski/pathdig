import glob
import os
import pandas as pd
import random

files = glob.glob("*neighborhood.tab")
outcols = ['sid']
dfs = list()
for f in files:
    print(f)
    structid,chain,_ = os.path.basename(f).split("_",2)
    cdf = pd.read_csv(f,sep="\t")
    if len(dfs)==0:
        outcols += [x for x in list(cdf.columns) if x!='variant']
    cdf['structid']=structid
    cdf['chain']=chain
    dfs.append(cdf)

final_df = pd.concat(dfs)
final_df['sid'] = final_df.structid+"_"+final_df.chain+"_"+final_df.variant

normcols = ['sid']+[x for x in outcols if "norm" in x or "delta" in x]
#normout = ['sid']+[x for x in outcols if "norm" in x or "total" in x or "delta" in x] #Include total counts?
final_df.to_csv("spatial_features.tab",sep="\t",columns=outcols,header=True,index=False)
final_df.to_csv("spatial_features_normonly.tab",sep="\t",columns=normcols,header=True,index=False)

# Randomly select 3 of each label from each structure
ixs = list()
for x in final_df.structid.unique():
    gnomad_ix = list(final_df[(final_df.structid==x) & (final_df.label=='gnomad')].index)
    gr = random.sample(gnomad_ix,min(len(gnomad_ix),3))
    clinvar_ix = list(final_df[(final_df.structid==x) & (final_df.label=='clinvar')].index)
    cr = random.sample(clinvar_ix,min(len(clinvar_ix),3))
    ixs += gr
    ixs += cr

final_df.iloc[ixs].to_csv("spatial_features_subset.tab",sep="\t",columns=outcols,header=True,index=False)    
final_df.iloc[ixs].to_csv("spatial_features_subset_normonly.tab",sep="\t",columns=normcols,header=True,index=False)        
