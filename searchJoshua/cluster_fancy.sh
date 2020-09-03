#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:00:00
#SBATCH --job-name=aispVS
#SBATCH --output=slurmOut/slurm_%a.out

# INPUT
# $1 Index of the job to run
# $2 full file name (including directory) for the data file

umask 077

dataDir="$2"
index=$1
job=$SLURM_JOB_ID

module purge
module load matlab/R2018b

cat<<EOF | matlab -nodisplay -nosplash
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

cluster_fcn_fancy("$dataDir", job_id, $index);

EOF
