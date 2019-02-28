
# Prepare a slurmfile for pathprox calculations

# Create new folder structure with directory for each structure - Done
# Create a command file in there which contains one pathprox command for each variant - Done
# Each Variant should get its own output folder within the structure folder - Done
# Also each variant needs its own input files for benign and pathogenic variants and a slightly modified version of its fasta file with a changed header - Done

import argparse
import logging
import sys
import os
#import shutil
import pandas as pd

ROOT_DIR = "/dors/capra_lab/users/sliwosgr/pathdig/pathprox/"
MAIL = "gregory.r.sliwoski@vanderbilt.edu"
version = 2.0

parser = argparse.ArgumentParser(description='Prepare file structure, data and slurm script for pathprox runs')
parser.add_argument('ddGprep_out', help="Folderstructure from a ddGprep.py run to reuse pdb and fasta files and variant positions.")
parser.add_argument('variants', help="csv file in dataCollecterComplete.py run format to collect structure names and chain IDs.")
parser.add_argument('pathprox_out', help="Path where a new folder structure should be created or old files overwritten.")
#parser.add_argument('--fasta_dir', '-f', type=str, help="Specify custom directory containing fastas in '<pdbid>.fasta' format.")

args = parser.parse_args()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

if not os.path.exists(args.pathprox_out):
    os.makedirs(args.pathprox_out)

fh = logging.FileHandler(os.path.join(args.pathprox_out,'pathproxPrepper.log'), 'w')
fh.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)

formatter = logging.Formatter('[%(asctime)s]%(levelname)s:%(message)s', '%m/%d/%Y %I:%M:%S %p')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(ch)

logger.info("Started pathproxPrepper.py v{}".format(version))
logger.info("User info: ROOT_DIR={}, MAIL={}".format(ROOT_DIR,MAIL))

#if args.fasta_dir == "":
#    logger.info("No custom fasta directory specified.")
#else:
#    logger.info("Fasta files will be searched in '" + args.fasta_dir + "'.")

# unpDict first collects pdbids as keys and provides information about chain and uniprot id
unpDict = {}
unpCounter = 0
variants = pd.read_csv(args.variants)
map_df = variants[['structid','chain','unp']].drop_duplicates()
#variants = open(args.variants, 'r')
unpCounter = map_df.shape[0]

#for line in variants:
#    if 'structid' in line:
#        pass
#    else:
#        line = line[:-1].split(',')
#        pdbid = line[0]
#        chain = line[1]
#        uniprot = line[2]
#        if pdbid not in unpDict:
#            unpCounter += 1
#            unpDict[pdbid] = [chain, uniprot]

#variants.close()
logger.info("Found " + str(unpCounter) + " potential structures.")

# Check for usable structures and their variants, create folder structure, create all needed csv files and cmd files, list structures for slurm
structList = []
unpCounter = 0
variantCounter = 0

with open(os.path.join(args.pathprox_out,"sliwosgr.config"),"w") as outfile:
        outfile.write("[UserSpecific]\n[SlurmParametersAll]\nmail-user={}\nmail-type=FAIL\n[SlurmParametersPathProxCOSMIC]\n".format(MAIL))
        outfile.write("time = 1-0\nmem = 30GB\n[SlurmParametersPathProxClinvar]\ntime = 1-0\nmem = 30GB\n[SlurmParametersPathProx]\ntime = 1-0\nmem = 30GB")

for index,row in map_df.iterrows():
    structid = row.structid
    chain = row.chain
    unp = row.unp
    if os.path.exists(args.ddGprep_out + structid + "/variants.csv"):
        unpCounter += 1
        structList.append(structid)
        current_out = os.path.join(args.pathprox_out,structid)
        if not os.path.isdir(current_out):
            os.makedirs(current_out)
        
        # keys are already in the right shape for pathprox variant format, value just tells clinvar or gnomad
#        variantDict = {}
#        print(args.ddGprep_out)
#        print(structid)
#        print(os.path.join(args.ddGprep_out,structid,"/variants.csv"))
        variants = pd.read_csv(os.path.join(args.ddGprep_out,structid,"variants.csv"))
        variants['variantID'] = variants.ref_amino_acid + variants.clean_pos.map(str) + variants.alt_amino_acid
        variantCounter += len(variants.variantID.unique())

        # Because pathprox needs fasta header to match uniprot exactly
        newfasta = ">{}|{}|hihihihihih\n".format(structid,unp)
        with open(os.path.join(args.ddGprep_out,structid,"{}_{}.fasta".format(structid,chain))) as infile:
            for line in infile.readlines()[1:]:
                newfasta += line.strip()
        newfasta += "\n"
        with open(os.path.join(args.pathprox_out,structid,"{}.fasta".format(structid)),'w') as outfile:
            outfile.write(newfasta)
#        shutil.copy(os.path.join(args.ddGprep_out,structid,"{}_{}.fasta".format(structid,chain)),
#                 os.path.join(args.pathprox_out, structid,"{}.fasta".format(structid)))

#        with open(os.path.join(args.ddGprep_out,structid,"/variants.csv"), 'r') as variants:
#            
#            for line in variants:
#                if 'structid' in line:
#                    pass
#                else:
#                    variantCounter += 1
#                    line = line[:-1].split(',')
#                    variantID = line[10]+line[9]+line[11]
#                    variantDict[variantID] = line[8]
        
        # Prepare fasta file
#        if args.fasta_dir == "":
#            oldfasta = open(args.ddGprep_out + structid + "/" + structid + "_" + unpDict[structid][0] + ".fasta", 'r')
#            sequence = oldfasta.readlines()[1]
#            oldfasta.close()
#            fasta = open(args.pathprox_out + structid + "/" + structid + ".fasta", 'w')
#            fasta.write(">" + structid + "|" + unpDict[structid][1] + "|Nonesense to fill the header\n")
#            fasta.write(sequence)
#            fasta.close()
#        else:
#            copyfile(args.fasta_dir + structid + ".fasta", args.pathprox_out + structid + "/" + structid + ".fasta")

        
        # Prepare config file
#        config = open(args.pathprox_out + "eisenhp.config", 'w')
#        config.write("[UserSpecific]\n[SlurmParametersAll]\nmail-user=eisenhuth451@gmail.com\nmail-type=all\n[SlurmParametersPathProxCOSMIC]\n")
#        config.write("time = 1-0\nmem = 30GB\n[SlurmParametersPathProxClinvar]\ntime = 1-0\nmem = 30GB\n[SlurmParametersPathProx]\ntime = 1-0\nmem = 30GB")
#        config.close()
               
        listoutfile = os.path.join(args.pathprox_out,structid,"{}_{}.csv")
       
        for vi,vr in variants.iterrows():
            # Create pathogenic and benign variant lists
            current_id = vr.variantID
            current_var = os.path.join(args.pathprox_out,structid,"{}_{}.csv")
            variants[(variants.variantID!=current_id) & (variants.label=="clinvar")].to_csv(
                                                                                    listoutfile.format("clinvar",current_id),
                                                                                    columns=['variantID'],
                                                                                    index=False,
                                                                                    header=False)
            variants[(variants.variantID!=current_id) & (variants.label=="gnomad")].to_csv(
                                                                                    listoutfile.format("gnomad",current_id),
                                                                                    columns=['variantID'],
                                                                                    index=False,
                                                                                    header=False)
            if variants[~variants.label.isin(["clinvar","gnomad"])].shape[0]>0:
                logger.error("Found variant unable to categorize.")            
#            clinvars = open(args.pathprox_out + structid + "/clinvar_" + variant + ".csv", 'w')
#            gnomads = open(args.pathprox_out + structid + "/gnomad_" + variant + ".csv", 'w')
#            for diffVariant in variantDict:
#                if diffVariant == variant:
#                    pass
#                else:
#                    if variantDict[diffVariant] == 'gnomad':
#                        gnomads.write(diffVariant + "\n")
#                    elif variantDict[diffVariant] == 'clinvar':
#                        clinvars.write(diffVariant + "\n")
#                    else:

#            clinvars.close()
#            gnomads.close()
            # Add command line to cmd file
            pdbfile = os.path.join(args.ddGprep_out, structid,"{}_{}.pdb".format(structid,chain))
            fastafile = "{}.fasta".format(structid)
            configfile = os.path.join(ROOT_DIR, args.pathprox_out,"sliwosgr.config")
            
            with open(os.path.join(args.pathprox_out, structid, "runpathprox.cmd"), 'a+') as cmdFile:
                cmdFile.write("pathprox2.py ")
                cmdFile.write("-c /dors/capra_lab/users/psbadmin/config/global.config ")
                cmdFile.write("--fasta {} ".format(fastafile))
                cmdFile.write("-u {} ".format(configfile))
                cmdFile.write("{} ".format(pdbfile))
                cmdFile.write("NM_00000.1 ")
                cmdFile.write("{} ".format(current_id))
                cmdFile.write("--chain={} ".format(chain))
                cmdFile.write("--neutral=gnomad_{}.csv ".format(current_id))
                cmdFile.write("--pathogenic=clinvar_{}.csv ".format(current_id))
                cmdFile.write("--radius=D ")
                cmdFile.write("--sqlcache={} ".format(os.path.join(ROOT_DIR,current_out,"{}_sql".format(current_id))))
                cmdFile.write("--overwrite ")
                cmdFile.write("--outdir {} ".format(os.path.join(ROOT_DIR,current_out,"{}_out".format(current_id))))
                cmdFile.write("--uniquekey {}\n\n".format(current_id))

logger.info("Created files for {} variants and {} structures.".format(variantCounter,unpCounter))

with open(os.path.join(args.pathprox_out,"pathprox.slurm"), 'w') as slurm:
    slurm.write("#!/bin/sh\n")
    slurm.write("#SBATCH --mail-user={}\n".format(MAIL))
    slurm.write("#SBATCH --mail-type=FAIL\n")
    slurm.write("#SBATCH --job-name=pathprox\n")
    slurm.write("#SBATCH --ntasks=1\n")
    slurm.write("#SBATCH --time=48:00:00\n")
    slurm.write("#SBATCH --mem=30G\n")
    maxRuns = len(structList)
    slurm.write("#SBATCH --array=0-{}\n".format(maxRuns-1))
    slurm.write("#SBATCH --output={}\n\n\n".format(os.path.join(ROOT_DIR,args.pathprox_out,"pathprox_%A_%a.out")))
    # export settings
    slurm.write("source /dors/capra_lab/users/psbadmin/psb_prep.bash\n\n")
    # Report outputs
    slurm.write("echo \"SLURM_JOBID: \" $SLURM_JOBID\n")
    slurm.write("echo \"SLURM_ARRAY_TASK_ID: \" $SLURM_ARRAY_TASK_ID\n")
    slurm.write("echo \"SLURM_ARRAY_JOB_ID: \" $SLURM_ARRAY_JOB_ID\n\n")
    # Make an array with all left to be used structures
    slurm.write("structures=({})\n".format(" ".join(structList)))
    slurm.write("activeStruct=${structures[${SLURM_ARRAY_TASK_ID}]}\n")
    slurm.write("echo \"ACTIVE_STRUCTURE: \" $activeStruct\n\n\n")
    # Try to get into the structure folder, exit if it is not possible
    sloppydir = os.path.join(ROOT_DIR,args.pathprox_out)
    if sloppydir[-1]!="/":
        sloppydir += "/"
    slurm.write("cd {}${{activeStruct}}/\n".format(sloppydir))
    slurm.write("if [ $? != 0 ]; then\n")
    slurm.write("echo \"Failure at script launch: Unable to change to directory\"\n")
    slurm.write("exit 1\n")
    slurm.write("fi\n\n\n")
    # Actually start pathprox
    slurm.write("source {}${{activeStruct}}/runpathprox.cmd > {}${{activeStruct}}/runpathprox.log\n\n\n".format(sloppydir,sloppydir))
    slurm.write("echo \"Finished program run.\"\n\n\n")

logger.info("Finished program run.")
