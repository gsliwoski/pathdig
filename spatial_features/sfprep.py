import sys
import os
import pandas as pd

try:
    variants = sys.argv[1]
    structures = sys.argv[2]

except:
    sys.exit("python sfprep.py variants.csv structure_path_root_dir")

assert os.path.isdir(structures), "{} not path".format(structures)
assert os.path.isfile(variants), "{} not found".format(variants)

vdf = pd.read_csv(variants)
vdf = vdf[['structid','chain','clean_pos','ref_amino_acid','alt_amino_acid','label']]
vdf['variantID'] = vdf.ref_amino_acid + map(str,vdf.clean_pos) + vdf.alt_amino_acid

for x in vdf.groupby(["structid","chain"]):
    s,varlines = x
    structid,chain = s
    pdbfile = os.path.join(structures,structid,"{}_{}.pdb".format(structid,chain))
    assert os.path.isfile(pdbfile), "pdbfile {} not found".format(pdbfile)
    varout = "{}_vars.ls".format(structid)
    varlines.to_csv(varout,header=False,index=False,sep=" ",columns=['variantID','label'])
    print "python /dors/capra_lab/users/sliwosgr/ml_tools/custom_features/neighborhood_features/gen_neighbor_features.py {} {} {}".format(pdbfile,chain,varout)  
    
