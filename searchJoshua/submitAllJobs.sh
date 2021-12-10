#!/bin/bash

# Submits all jobs which don't have an associated results file.

# INPUT
# $1 Which jobs should the function run (if they have not already been 
#   run)? This number should be the begining of the range.
# $2 Which jobs should the function run (if they have not already been 
#   run)? This number should be the end of the range.
# $3 Current cluster ("main" or "nyu")

set -e

startJobs=$1
endJobs=$2
currentCluster=$3

# Check we are in the "searchJoshua" directory
currentDir=$(basename $(pwd))
if [[ $currentDir != searchJoshua ]]; then
    echo "Must be run from the subdirectory searchJoshua"
    exit 1
fi

dataDir="$(pwd)/Data/StandardFormat_participantExcluded.mat"

if [[ $currentCluster == main ]]; then
    launchScript="./cluster_fancy.sh"
elif [[ $currentCluster == nyu ]]; then
    launchScript="./cluster_fancy_nyu.sh"
else
    echo "Unrecognised cluster selected"
    exit 1
fi

# Loop through all job files in the direcotry
for idx in $(seq $startJobs $endJobs); do

    resultExists=0
    resultToFind=/pars/pars_$idx

    for filename in ./pars/*.mat; do
        if [[ "$filename" == *"$resultToFind"_* ]]; then
            resultExists=1
        fi
    done
    
    if [[ $resultExists == 0 ]]; then
        sbatch $launchScript "$idx" "$dataDir"
    else
        echo "Results for job index $idx found, so skipping job submission"
    fi
done