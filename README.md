Paul's datasets:

1. ddG monomer 
* I did not rerun these calculations but copied them over into ddGprep_out0. The only thing that is required from here are the pdb structures used in other scripts. = /dors/capra_lab/users/sliwosgr/pathdig/ddGprep_out0
2. pathprox
* I identified potentially what was causing issues in Paul's original pathprox analysis (see below)
* I refactored his pathprox scripts and so to set up all the scripts and SLURM files, I run

python3 pathproxPrepper_gs.py /dors/capra_lab/users/sliwosgr/pathdig/ddGprep_out0/ /dors/capra_lab/users/sliwosgr/pathdig/variants/resultUnique.csv PP_run

* This generates a directory called PP_run containing all necessary files.
* I then submit PP_run/pathprox.slurm as a batch array with the number of structures in the array
* Then run `pathproxCollector.py ddGscores.csv PP_run` to generate a summary of the pathprox runs once they complete.
* ddGscores.csv is the output from Paul's ddG collector and the pathprox scores are attached to these rows per variant.
* By default, this creates a file called scoredVariants.csv whre two additional columns have been added: 'variantID' and 'pathprox'

3. spatial neighborhood descriptor
* Located in spatial_descriptors/
* run sfprep.py with the csv of variants and the path to the ddG structures (ddGprep_out0) full path.
* This generates all of the necessary files
* gen_features.cmd contains all the lines to generate the features across variants. This does not take long so it doesn't need to be submitted.

Example line:

`python /dors/capra_lab/users/sliwosgr/ml_tools/custom_features/neighborhood_features/gen_neighbor_features.py /dors/capra_lab/users/sliwosgr/pathdig/ddGprep_out0/5FQD/5FQD_B.pdb B 5FQD_vars.ls`

* Once that is finished, run collect_sf.py to gather together all the spatial features across the tab files. These are the spatial features used to train models. They can be attached to ddG and pathprox features or any other feature set that is tab delimited and must contain a column of the variant id (by default, just use sid for the column name)
* In random_forest/ there is some files for running random forest predictions. By default everything goes into output/ and the files should be pretty self explanatory for all the different metrics of a rf classifier.
* Any set of features can be used here as long as they are formatted as above. For the labels, just use sf_y.tab for all the runs unless you end up generating your own set of variants in which case you need a labels file that looks like sf_y.tab but with your own id's and labels.
* run `/dors/capra_lab/users/sliwosgr/ml_tools/random_forest.py -f sf_x.tab -l sf_y.tab` for the example run. It's pretty fast so currently don't need to submit to accre. There are a bunch of other flags and ability to set your own parameters with a params file. Right now, the default parameters have not undergone any hyperoptimization and are likely not the best to use, but are what I've been using to do initial test runs and feature creation.

* You can also run `/dors/capra_lab/users/sliwosgr/ml_tools/ann.py -f sf_x.tab -l sf_y.tab` to do the same thing using a 2 hidden layer ANN instead of random forest. This can also take parameters and the default ones were just arbitrarily set and have not undergone any hyperoptimization so defining params will likely get better performance.



Problems with Pathprox:

Previously:

What I currently think happened:

1. Paul defined variants in terms to seqid position (as in the position on the transcript/uniprot)
2. Paul defines fasta sequence initially in terms of cleaned PDB fasta sequence
2ex: structid = 1A4Y, chain = B, variant = C39W
/dors/capra_lab/users/eisenhp/pathproxPrep/1A4Y/C39W_out
In this case, the seqid = 39 and the clean_pos = 39
3. Paul grabs the raw PDB file.
4. Paul has problems because the cleaned pdb fasta may be missing residues that were removed because they were missing backbone atoms or any other potential changes to the number of residues in the structure. It may have been that one of the variants fell on a residue in the raw PDB fasta that was lost during the cleaning to the cleaned fasta. This would cause errors when attempting to convert numbering.
5. Chris provides Paul with reference sequences that contain the entire CURRENT UNIPROT sequences.
-This is the critical part.
In the example for 1A4Y, the uniprot is P03950.
In the UNIPROT sequence, Asp is at residue 39 and Cys is actually at 63. The first 24 residues of the Uniprot sequence are missing from the PDB sequence. And since the seqid column and clean_pos column are both 39, I'm assuming the transcript aligned to this PDB was also missing these 24 residues. In this case, the numbering between Uniprot and initial seqid, which I assume comes from PDBMap, is NOT the same as the uniprot.
6. pathprox aligns the raw PDB (Cys at 37) with the Uniprot sequence (Cys at 64).
7. pathprox assumes that Paul gave the variant 37 as the reference position, which should have been 64.
8. pathprox correctly translates 37 into original PDB numbering of 13 based on the assumption that Paul is interested in position 37 of the uniprot sequence.
-This is the second critical part
pathprox does not check to ensure that your reference AA matches the AA found on the reference sequence. In this case, C!=D which should have been flagged or a strong warning issued which would have revealed that the fastas provided were problematic.

This probably didn't happen in all structures and looking at a few examples it can get very confusing especially since it took a while to understand what pathprox was doing in terms of renumbering things.

To fix this, just make sure there is consistency in fasta vs PDB file vs variant numbering.

Since the cleaned PDB was given to ddG monomer, it's best to use the cleaned PDB files in addition to the cleaned PDB fasta and use the clean_pos column for all variant position numbering.


Additionally, it should be noted that Paul's cutoff was 3 of each type of variant for structure selection.

This means that with loo pathprox, it can not run with any structures where by leaving one of three out results in 2 remaining for that condition.

Other observations:

1FYH may be an example where ddG provides info that pathprox lacks.
V33G is clinvar while V33I is gnomad. These are at the same position which would give them the same PP score. However, V33G has a ddG of 8.2 while V33I has ddG of 1.503.
However, this also points out a potential problem, their PP scores in this training set are not actually the same since they are based on a leave-one-out method which changes the number of clinvar and gnomad variants at that position depending on which is left out.

Normalized features should likely be used to avoid the different number of variants in the leave-one-out training strategy that fall within the same structure.
