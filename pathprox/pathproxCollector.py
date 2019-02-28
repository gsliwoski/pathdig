# Read the ddG score csv
# Gather all the pathprox scores
# Merge with ddG scores
# write out merged file and report any missing scores

import argparse
import pandas as pd
import os

import sys
parser = argparse.ArgumentParser(description='Collect pathprox scores and append to ddG scores')
parser.add_argument('ddGs', help="Inputfile containing all ddG scores from a ddGcollector.py run.")
parser.add_argument('pathprox', help="Folderstructure from a pathprox run prepared by pathproxPrepper.py")
parser.add_argument('--output', '-o', type=str, default='scoredVariants.csv', help='File to write out results.')

args = parser.parse_args()

ddGs = pd.read_csv(args.ddGs)
ddGs['variantID'] = ddGs.ref_amino_acid + ddGs.clean_pos.map(str) + ddGs.alt_amino_acid
print("Retrieving {} variant pathprox scores".format(ddGs.shape[0]))
if "pathprox" in ddGs.columns:
    ddGs.drop("pathprox",axis=1,inplace=True)
outcols = list(ddGs.columns) + ["pathprox"]
all_scores = None
for sc in ddGs[['structid','chain']].drop_duplicates().iterrows():
    structid = sc[1].structid
    chain = sc[1].chain
    print("collecting {}_{}".format(structid,chain))
    for var in ddGs[ddGs.structid==structid].variantID:
        scorefile = os.path.join(args.pathprox, structid, "{}_out".format(var),"{}_{}_{}_D_summary.csv".format(structid,chain,chain))
        if not os.path.isfile(scorefile):
            print("{} not found".format(scorefile))
            continue
        try:
            current_score = pd.read_csv(scorefile,sep="\t")
        except:
            continue            
        mapper = {"{}_pathprox".format(var):"pathprox",
                  "{}_pathcon".format(var):"pathcon",
                  "{}_neutcon".format(var):"neutcon"}
        current_score.rename(mapper,inplace=True,axis=1)        
        if current_score[~current_score.pathprox.isnull()].shape[0] == 0:
            continue
        current_score['structid'] = structid
        current_score['chain'] = chain
        current_score['variantID'] = var            
        if all_scores is None:
            all_scores = current_score
        else:
            all_scores = pd.concat([all_scores, current_score])
if all_scores is None or all_scores.shape[0] == 0:
    sys.exit("failed to get any pathprox scores")
else:
    print("Retrived {} pathprox scores".format(all_scores.shape[0]))

ddGs = ddGs.merge(all_scores,on=["structid","chain","variantID"],how="left")
missingpp = ddGs[ddGs.pathprox.isnull()]
print("{} structures and {} variants have missing pathprox scores, see missing_pp.csv".format(len(missingpp.structid.unique()),
                                                                                              missingpp.shape[0]))
missingpp.sort_values(["structid","chain","variantID"]).to_csv("missing_pp.csv",index=False,header=True,columns=outcols[:-1])

print("Writing final scorefile to {}".format(args.output))
ddGs.sort_values(["structid","chain","variantID"]).to_csv(args.output,index=False,header=True,columns=outcols)
