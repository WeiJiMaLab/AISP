#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=7GB
#SBATCH --job-name=refit_idx
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=re_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/bigAISP/detectionAspen'))

% fixed model fitting settings
load('idxstocomplete.mat')
runs = csvread('runs.csv');

idxx = $SLURM_ARRAY_TASK_ID;
idx = idxlist(idxx);
imodel = runs(idx, 3);
isubj = runs(idx, 1);
irep = runs(idx, 2);

fit_cluster_ibs(imodel,isubj,irep)

EOF
