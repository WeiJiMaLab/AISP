#!/bin/bash
#SBATCH --job-name=aispVS
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:00:00
#SBATCH --output=slurmOut/slurm_%j.out
#SBATCH --export=NONE

# INPUT
# $1 Index of the job to run
# $2 full file name (including directory) for the data file

umask 077
set -e

source /sw/batch/init.sh

dataDir="$2"
index=$1
job=$SLURM_JOB_ID

module load matlab/2018b

cat<<EOF | matlab -nodisplay -nosplash
job_id = str2num(strjoin(regexp('$job','\d','match'), ''))
rng(job_id)

tmpFolder = fullfile("$RRZ_LOCAL_TMPDIR", 'clusterJobs', ...
    num2str(job_id), "$index")
disp('Folder in use for temporary storage of files...')
disp(tmpFolder)
disp('')
mkdir(tmpFolder)

thisClst = parcluster();
thisClst.JobStorageLocation = tmpFolder;
parpool(thisClst, [1, 32])

# cluster_fcn_fancy("$dataDir", job_id, $index);

EOF
