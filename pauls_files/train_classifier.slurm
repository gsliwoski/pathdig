#!/bin/sh
#SBATCH --mail-user=eisenhuth451@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=classifier
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --array=0-1
#SBATCH --output=/dors/capra_lab/users/eisenhp/classifiers/logfile_classifier_%A_%a.out
#SBATCH --partition=debug

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

cd /dors/capra_lab/users/eisenhp/classifiers/
if [ $? != 0 ]; then
echo "Failure at script launch: Unable to change to directory"
exit 1
fi

python3 classificationTrainer.py --max_iter 10000 --save_interval 1000 --report_interval 500 --activation relu annotatedScores.csv $SLURM_ARRAY_TASK_ID

