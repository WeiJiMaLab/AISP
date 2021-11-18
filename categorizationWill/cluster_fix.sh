#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=fixWill
#SBATCH --mail-type=END
#SBATCH --mail-user=hhs4@nyu.edu
#SBATCH --output=slurm-will/slurm_fix_%j.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID
ppn=$SLURM_JOB_CPUS_PER_NODE
module purge
module load matlab/2020b
export MATLABPATH=$HOME/matlab-output

cat<<EOF | matlab -nodisplay
cd ~/AISP/categorizationWill
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

tmpfolder = sprintf('/state/partition1/job-%d/.matlab/', job_id)
mkdir(tmpfolder)
clust = parcluster();
clust.JobStorageLocation = tmpfolder;
parpool('threads')

cluster_fcn_fix(job_id);

EOF
