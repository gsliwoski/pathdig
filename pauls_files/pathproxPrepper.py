
# Prepare a slurmfile for pathprox calculations

# Create new folder structure with directory for each structure - Done
# Create a command file in there which contains one pathprox command for each variant - Done
# Each Variant should get its own output folder within the structure folder - Done
# Also each variant needs its own input files for benign and pathogenic variants and a slightly modified version of its fasta file with a changed header - Done

import argparse
import logging
import sys
import os
from shutil import copyfile

version = 1.1

parser = argparse.ArgumentParser(description='Prepare file structure, data and slurm script for pathprox runs')
parser.add_argument('ddGprep_out', help="Folderstructure from a ddGprep.py run to reuse pdb and fasta files and variant positions.")
parser.add_argument('variants', help="csv file in dataCollecterComplete.py run format to collect structure names and chain IDs.")
parser.add_argument('pathprox_out', help="Path where a new folder structure should be created or old files overwritten.")
parser.add_argument('--fasta_dir', '-f', type=str, help="Specify custom directory containing fastas in '<pdbid>.fasta' format.")

args = parser.parse_args()


logger = logging.getLogger()
logger.setLevel(logging.INFO)

if not os.path.exists(args.pathprox_out):
    os.makedirs(args.pathprox_out)

fh = logging.FileHandler(args.pathprox_out + 'pathproxPrepper.log', 'w')
fh.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)

formatter = logging.Formatter('[%(asctime)s]%(levelname)s:%(message)s', '%m/%d/%Y %I:%M:%S %p')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(ch)


logger.info("Started pathproxPrepper.py v" + str(version))

if args.fasta_dir == "":
    logger.info("No custom fasta directory specified.")
else:
    logger.info("Fasta files will be searched in '" + args.fasta_dir + "'.")

# unpDict first collects pdbids as keys and provides information about chain and uniprot id
unpDict = {}
unpCounter = 0
variants = open(args.variants, 'r')
for line in variants:
    if 'structid' in line:
        pass
    else:
        line = line[:-1].split(',')
        pdbid = line[0]
        chain = line[1]
        uniprot = line[2]
        if pdbid not in unpDict:
            unpCounter += 1
            unpDict[pdbid] = [chain, uniprot]

variants.close()
logger.info("Found " + str(unpCounter) + " potential structures.")

# Check for usable structures and their variants, create folder structure, create all needed csv files and cmd files, list structures for slurm
structList = []
unpCounter = 0
variantCounter = 0
for structid in unpDict:
    if os.path.exists(args.ddGprep_out + structid + "/variants.csv"):
        unpCounter += 1
        structList.append(structid)
        if not os.path.exists(args.pathprox_out + structid + "/"):
            os.makedirs(args.pathprox_out + structid + "/")
        
        # keys are already in the right shape for pathprox variant format, value just tells clinvar or gnomad
        variantDict = {}
        variants = open(args.ddGprep_out + structid + "/variants.csv", 'r')
        for line in variants:
            if 'structid' in line:
                pass
            else:
                variantCounter += 1
                line = line[:-1].split(',')
                variantID = line[10]+line[9]+line[11]
                variantDict[variantID] = line[8]
        variants.close()
        
        # Prepare fasta file
        if args.fasta_dir == "":
            oldfasta = open(args.ddGprep_out + structid + "/" + structid + "_" + unpDict[structid][0] + ".fasta", 'r')
            sequence = oldfasta.readlines()[1]
            oldfasta.close()
            fasta = open(args.pathprox_out + structid + "/" + structid + ".fasta", 'w')
            fasta.write(">" + structid + "|" + unpDict[structid][1] + "|Nonesense to fill the header\n")
            fasta.write(sequence)
            fasta.close()
        else:
            copyfile(args.fasta_dir + structid + ".fasta", args.pathprox_out + structid + "/" + structid + ".fasta")
        
        # Prepare config file
        config = open(args.pathprox_out + "eisenhp.config", 'w')
        config.write("[UserSpecific]\n[SlurmParametersAll]\nmail-user=eisenhuth451@gmail.com\nmail-type=all\n[SlurmParametersPathProxCOSMIC]\n")
        config.write("time = 1-0\nmem = 30GB\n[SlurmParametersPathProxClinvar]\ntime = 1-0\nmem = 30GB\n[SlurmParametersPathProx]\ntime = 1-0\nmem = 30GB")
        config.close()
        
        cmdFile = open(args.pathprox_out + structid + "/runpathprox.cmd", 'w')
        
        for variant in variantDict:
            # Create pathogenic and benign variant lists
            clinvars = open(args.pathprox_out + structid + "/clinvar_" + variant + ".csv", 'w')
            gnomads = open(args.pathprox_out + structid + "/gnomad_" + variant + ".csv", 'w')
            for diffVariant in variantDict:
                if diffVariant == variant:
                    pass
                else:
                    if variantDict[diffVariant] == 'gnomad':
                        gnomads.write(diffVariant + "\n")
                    elif variantDict[diffVariant] == 'clinvar':
                        clinvars.write(diffVariant + "\n")
                    else:
                        logger.error("Found variant unable to categorize.")
            clinvars.close()
            gnomads.close()
            # Add command line to cmd file
            cmdFile.write("pathprox2.py ")
            cmdFile.write("-c /dors/capra_lab/users/psbadmin/config/global.config ")
            cmdFile.write("--fasta " + structid + ".fasta ")
            cmdFile.write("-u /dors/capra_lab/users/eisenhp/" + args.pathprox_out + "eisenhp.config ")
            cmdFile.write("/dors/capra_lab/users/eisenhp/" + args.ddGprep_out + structid + "/" + structid + ".pdb ")
            cmdFile.write("NM_00000.1 ")
            cmdFile.write(variant + " ")
            cmdFile.write("--chain=" + unpDict[structid][0] + " ")
            cmdFile.write("--neutral=gnomad_" + variant + ".csv ")
            cmdFile.write("--pathogenic=clinvar_" + variant + ".csv ")
            cmdFile.write("--radius=D ")
            cmdFile.write("--sqlcache=/dors/capra_lab/users/eisenhp/" + args.pathprox_out + structid + "/" + variant + "_sql/ ")
            cmdFile.write("--overwrite ")
            cmdFile.write("--outdir /dors/capra_lab/users/eisenhp/" + args.pathprox_out + structid + "/" + variant + "_out/ ")
            cmdFile.write("--uniquekey " + variant + "\n\n")
        
        cmdFile.close()

logger.info("Created files for " + str(variantCounter) + " variants and " + str(unpCounter) + " structures.")

slurm = open(args.pathprox_out + "pathprox.slurm", 'w')
slurm.write("#!/bin/sh\n")
slurm.write("#SBATCH --mail-user=eisenhuth451@gmail.com\n")
slurm.write("#SBATCH --mail-type=ALL\n")
slurm.write("#SBATCH --job-name=pathprox\n")
slurm.write("#SBATCH --ntasks=1\n")
slurm.write("#SBATCH --time=48:00:00\n")
slurm.write("#SBATCH --mem=30G\n")
maxRuns = len(structList)
slurm.write("#SBATCH --array=0-" + str(maxRuns-1) + "\n")
slurm.write("#SBATCH --output=/dors/capra_lab/users/eisenhp/" + args.pathprox_out + "pathprox_%A_%a.out\n\n\n")
# export settings
slurm.write("source /dors/capra_lab/users/psbadmin/psb_prep.bash\n\n")
# Report outputs
slurm.write("echo \"SLURM_JOBID: \" $SLURM_JOBID\n")
slurm.write("echo \"SLURM_ARRAY_TASK_ID: \" $SLURM_ARRAY_TASK_ID\n")
slurm.write("echo \"SLURM_ARRAY_JOB_ID: \" $SLURM_ARRAY_JOB_ID\n\n")
# Make an array with all left to be used structures
slurm.write("structures=(" + " ".join(structList) + ")\n")
slurm.write("activeStruct=${structures[${SLURM_ARRAY_TASK_ID}]}\n")
slurm.write("echo \"ACTIVE_STRUCTURE: \" $activeStruct\n\n\n")
# Try to get into the structure folder, exit if it is not possible
slurm.write("cd /dors/capra_lab/users/eisenhp/" + args.pathprox_out + "${activeStruct}/\n")
slurm.write("if [ $? != 0 ]; then\n")
slurm.write("echo \"Failure at script launch: Unable to change to directory\"\n")
slurm.write("exit 1\n")
slurm.write("fi\n\n\n")
# Actually start ddG_monomer
slurm.write("source /dors/capra_lab/users/eisenhp/" + args.pathprox_out + "${activeStruct}/runpathprox.cmd > /dors/capra_lab/users/eisenhp/" + args.pathprox_out + "${activeStruct}/runpathprox.log\n\n\n")
slurm.write("echo \"Finished program run.\"\n\n\n")

slurm.close()

logger.info("Finished program run.")


