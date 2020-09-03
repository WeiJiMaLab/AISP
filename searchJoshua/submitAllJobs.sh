#!/bin/bash

# Submits all jobs which don't have an associated results file.

# INPUT
# $1 How many jobs should the function run (if they have not already been run)?
# $2 full file name (including directory) for the data file

numJobs=$1
dataDir=$2

# Loop through all job files in the direcotry
for idx in $(seq 1 $numJobs); do

    resultExists=0
    resultToFind=/pars/pars_$idx

    for filename in ./pars/*.mat; do
        if [[ "$filename" == *"$resultToFind"* ]]; then
            resultExists=1
        fi
    done
    
    if [[ $resultExists == 0 ]]; then
        sbatch ./cluster_fancy.sh "$idx" "$dataDir"
    else
        echo "Job index $idx found"
    fi
done