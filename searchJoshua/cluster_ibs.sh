#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=fitYanli
#SBATCH --mail-type=END
#SBATCH --mail-user=hhs4@nyu.edu
#SBATCH --output=slurm-Yanli/slurm_%a.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID
ppn=$SLURM_JOB_CPUS_PER_NODE
module purge
module load matlab/2018b
export MATLABPATH=$HOME/matlab-output

cat<<EOF | matlab -nodisplay
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

newdir = '$SCRATCH/cluster$job';

mkdir(newdir);
cluster_fcn_ibs(job_id,$index);

rmdir(newdir,'s')


EOF
