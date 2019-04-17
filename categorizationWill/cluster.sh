#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:02:00
#SBATCH --mem=2GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=hhs4@nyu.edu
#SBATCH --output=slurm_%j.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID
ppn=$SLURM_JOB_CPUS_PER_NODE
module purge
module load matlab/2018b
export MATLABPATH=$HOME/matlab-output

cat<<EOF | matlab -nodisplay
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)
addpath(genpath('$HOME/AISP'))

newdir = '$SCRATCH/cluster';

mkdir(newdir);
cluster{$index} = parcluster('local');
cluster{$index}.JobStorageLocation = newdir;
parpool(cluster{$index}, $ppn);
cluster_fcn(job_id);

rmdir(newdir,'s')


EOF