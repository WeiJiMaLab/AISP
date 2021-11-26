#!/bin/bash

# Submits all jobs which don't have an associated results file.

# INPUT
# $1 Which jobs should the function run (if they have not already been 
#   run)? This number should be the begining of the range.
# $2 Which jobs should the function run (if they have not already been 
#   run)? This number should be the end of the range.
# $3 full file name (including directory) for the data file

set -e

startJobs=$1
endJobs=$2
dataDir=$3

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
        sbatch ./cluster_fancy.sh "$idx" "$dataDir"
    else
        echo "Job index $idx found"
    fi
done