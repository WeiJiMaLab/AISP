#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --job-name=aispVS
#SBATCH --output=slurmOut/slurm_%a.out

# INPUT
# $1 full file name (including directory) for the data file

umask 077

dataDir="$1"
index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID

module purge
module load matlab/R2018b

cat<<EOF | matlab -nodisplay -nosplash
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

cluster_fcn_ibs("$dataDir", job_id, $index);

EOF
