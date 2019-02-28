/dors/capra_lab/users/sliwosgr/pathdig

(Need to untarzip the ddGprep_out0 before running pathprox prepper script)

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

Missing without good reason:
3QIC: E273Q

Missing with reason:
1A4Y: only 3 clinvar - no clinvar pathprox
1ALD: only 3 clinvar - no clinvar pathprox
1ALY: only 3 gnomad - no gnomad pathprox
1BZ1: only 3 gnomad- no gnomad pathprox
1DXT: only 3 gnomad - no gnomad pathprox
1FGG: only 3 clinvar - no clinvar pathprox
1GR3: only 3 gnomad - no gnomad pathprox
1HZD: only 3 gnomad - no gnomad pathprox
1IK9: only 3 clinvar - no clinvar pathprox
1L4D: only 3 clinvar - no clinvar pathprox
1LMT: only 3 gnomad - no gnomad pathprox
1LYW: only 3 gnomad - no gnomad pathprox
1NNL: only 3 clinvar - no clinvar pathprox
1NST: only 3 gnomad - no gnomad pathprox
1NTG: only 3 clinvar - no clinvar pathprox
1OIZ: only 3 gnomad - no gnomad pathprox
1PV8: only 3 clinvar - no clinvar pathprox
1SI4: only 3 gnomad - no gnomad pathprox
1T94: only 3 clinvar - no clinvar pathprox
1TJJ: only 3 clinvar - no clinvar pathprox
1Z6Z: only 3 gnomad - no gnomad pathprox
2EYI: only 3 gnomad - no gnomad pathprox
2FIM: only 3 gnomad - no gnomad pathprox
2OBV: only 3 gnomad - no gnomad pathprox
2ODU: only 3 clinvar - no clinvar pathprox
2VKQ: only 3 clinvar - no clinvar pathprox
2VT5: only 3 gnomad - no gnomad pathprox
2WUL: only 3 of either - no variants pathprox
2XT3: only 3 clinvar - no clinvar pathprox
2XZE: only 3 of either - no variants pathprox
2ZEJ: only 3 gnomad - no gnomad pathprox
3A1F: only 3 gnomad - no gnomad pathprox
3BFN: only 3 of either - no variants pathprox
3BH7: only 3 clinvar - no clinvar pathprox
3E77: only 3 clinvar - no clinvar pathprox
3EXE: only 3 gnomad - no gnomad pathprox
3FS1: only 3 gnomad - no gnomad pathprox
3G2F: only 3 gnomad - no gnomad pathprox
3H3G: only 3 gnomad - no gnomad pathprox
3KJP: only 3 gnomad - no gnomad pathprox
3KWB: only 3 of either - no variants pathprox
3L81: only 3 clinvar - no clinvar pathprox
3LC3: only 3 gnomad - no gnomad pathprox
3LQ3: only 3 clinvar - no clinvar pathprox
3LWK: only 3 clinvar - no clinvar pathprox
3QK3: only 3 clinvar - no clinvar pathprox
3R8Q: only 3 clinvar - no clinvar pathprox
3RFE: only 3 of either - no variants pathprox
3RMU: only 3 clinvar - no clinvar pathprox
3ZQS: only 3 of either - no variants pathprox
3ZYC: only 3 clinvar - no clinvar pathprox
4A9C: only 3 gnomad - no gnomad pathprox
4APO: only 3 gnomad - no gnomad pathprox
4AY9: only 3 gnomad - no gnomad pathprox
4B7L: only 3 gnomad - no gnomad pathprox
4BSJ: only 3 clinvar - no clinvar pathprox
4CC0: only 3 gnomad - no gnomad pathprox
4CH8: only 3 gnomad - no gnomad pathprox
4D0E: only 3 gnomad - no gnomad pathprox
4DM9: only 3 of either - no variants pathprox
4FVL: only 3 clinvar - no clinvar pathprox
4FVQ: only 3 gnomad - no gnomad pathprox
4GS7: only 3 gnomad - no gnomad pathprox
4GT4: only 3 clinvar - no clinvar pathprox
4HL4: only 3 of either - no variants pathprox
4IDQ: only 3 gnomad - no gnomad pathprox
4JA8: only 3 of either - no variants pathprox
4P5W: only 3 gnomad - no gnomad pathprox
4QAM: only 3 clinvar - no clinvar pathprox
4TTB: only 3 clinvar - no clinvar pathprox
4WND: only 3 gnomad - no gnomad pathprox
4WNU: only 2 clinvar - no clinvar pathprox
4XO7: only 3 clinvar - no clinvar pathprox
5A0C: only 3 gnomad - no gnomad pathprox
5BON: only 3 of either - no variants pathprox
5E26: only 3 gnomad - no gnomad pathprox
5KWY: only 3 gnomad - no gnomad pathprox
5LGD: only 3 clinvar - no clinvar pathprox

Missing from Paul's only:
1APZ: K13E
1BZ4: L7P, R4C
1FYH: C49Y, G141R, I59T, T3I, V33G
1FZC: A185T, G250D, L203R, L22Q, R16C (these make up the full clinvar variant set)
1G73: All variants
1HGU: A12T, A12V, F24Y, R15C, R18H
1HP7: L23P, R21C, R21H, S18R
1HWG: E13K, R12G
1IVY: Q21R
1L4D: G18R
1LYW: All variants
1NTG: All variants
1OOK: V15A
1SI4: V1A
1TJJ: All variants
1WSR: A20V, G16W, H11R, R3H
2A01: P3R, P4R, R10L
2BFF: Q3H
2BV5: A139V, E10G, G241S, R4H, S165N, Y83C
2FIM: All variants
2H2N: V61L, W41C, Y7H
2ODU: All variants
2PUQ: V6M
2QTZ: I60T, S23F, T56M
2R0N: R42C
3A1F: All variants
3BPT: All variants
3DKB: All variants
3ECR: A38S, H220N
3F7B: N150K, R192S, R212Q, R50H, V207I
3FJO: All variants
3IAR: All variants
3KCG: I5N, R22C
3LC3: All variants
3PDF: All variants
3R8Q: All variants
3TEG: S11C
3UP1: G19V, S20L, S20W
3ZE2: R39W
4ANP: Clinvars
4CG4: All variants
4CH8: M17T
4DO4: All variants
4GS7: G61D
4H4F: R8Q
4OKH: All variants
4P0J: All variants
4TTB: All variants
4X30: N10D
4X4W: All variants
5A5I: All variants
5BV7: F243L
5E2P: All variants
5KWY: All variants

Other observations:

1FYH may be an example where ddG provides info that pathprox lacks.
V33G is clinvar while V33I is gnomad. These are at the same position which would give them the same PP score. However, V33G has a ddG of 8.2 while V33I has ddG of 1.503.
However, this also points out a potential problem, their PP scores in this training set are not actually the same since they are based on a leave-one-out method which changes the number of clinvar and gnomad variants at that position depending on which is left out.
