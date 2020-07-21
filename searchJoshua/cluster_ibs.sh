#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --job-name=aispVS
#SBATCH --output=slurm/slurm_%a.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID

module purge
module load matlab/2018b

cat<<EOF | matlab -nodisplay
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

newdir = 'tmpData/cluster$job';
mkdir(newdir);

cluster_fcn_ibs(job_id,$index);

rmdir(newdir,'s')


EOF
