#!/bin/bash
#SBATCH --job-name=aispVS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=slurmOut/slurm_%j.out
#SBATCH --export=NONE

# INPUT
# $1 Index of the job to run
# $2 full file name (including directory) for the data file

umask 077
set -e

dataDir="$2"
index=$1
job=$SLURM_JOB_ID

module purge
module load matlab/2018b

TMP="$TMPDIR"
echo 'System temp dirs:'
echo "$TMPDIR"
echo "$TMP"

export MATLAB_PREFDIR=$TMPDIR/.matlab/R2018b/
mkdir $TMPDIR/.matlab
cp -r $HOME/.matlab/R2018b $TMPDIR/.matlab

cat<<EOF | matlab -nodisplay -nosplash
job_id = str2num(strjoin(regexp('$job','\d','match'), ''));
rng(job_id)

tmpFolder = fullfile("$TMPDIR", 'clusterJobs', ...
    num2str(job_id), "$index");
disp('Folder in use for temporary storage of files...')
disp(tmpFolder)
disp('')
mkdir(tmpFolder)

thisClst = parcluster();
thisClst.JobStorageLocation = tmpFolder;
parpool(thisClst, [1, 32])
cluster_fcn_fancy("$dataDir", job_id, $index);
EOF
