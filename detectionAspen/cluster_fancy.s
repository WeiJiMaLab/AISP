#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=clusterfancy
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=o_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/bigAISP/detectionAspen'))

ttime = tic;
idx = $SLURM_ARRAY_TASK_ID;
cluster_fcn_fancy(idx)

fprintf('job complete in %f minutes \n',toc(ttime)/60)

EOF
