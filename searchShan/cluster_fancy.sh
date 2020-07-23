#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=fitShan
#SBATCH --mail-type=END
#SBATCH --mail-user=hhs4@nyu.edu
#SBATCH --output=slurm-Shan/slurm_%a.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID
ppn=$SLURM_JOB_CPUS_PER_NODE
module purge
module load matlab/2018b
export MATLABPATH=$HOME/matlab-output

cat<<EOF | matlab -nodisplay
cd ~/AISP/searchShan
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)
cluster_fcn_fancy(job_id,$index);

EOF
