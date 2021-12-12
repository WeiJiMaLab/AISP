#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=36:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=fitAspen
#SBATCH --mail-type=END
#SBATCH --mail-user=hhs4@nyu.edu
#SBATCH --output=slurm-aspen/slurm_%a.out

index=$SLURM_ARRAY_TASK_ID
job=$SLURM_JOB_ID
ppn=$SLURM_JOB_CPUS_PER_NODE
module purge
module load matlab/2020b
export MATLABPATH=$HOME/matlab-output

export MATLAB_PREFDIR=$TMPDIR/.matlab/R2020b/
mkdir $TMPDIR/.matlab
cp -r $HOME/.matlab/R2020b $TMPDIR/.matlab
mkdir $TMPDIR/.matlab/local_cluster_jobs
mkdir $TMPDIR/.matlab/local_cluster_jobs/R2020b

cat<<EOF | matlab -nodisplay
cd ~/AISP/detectionAspen

tic;
tmpfolder = sprintf('$TMPDIR/.matlab/local_cluster_jobs/R2020b')
clust = parcluster();
clust.JobStorageLocation = tmpfolder;
parpool(clust)

cluster_fcn_fancy($index)

ttime = toc;
fprintf('job complete in %f minutes \n',toc(ttime)/60);

EOF
