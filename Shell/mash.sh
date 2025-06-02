#!/bin/bash

#SBATCH --job-name=mash_screen
#SBATCH --account=xx
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

## Variables
DB=/PATH/plsdb_sketch.msh

echo 'Usage: sbatch mash.sh'

## Load modules
echo 'Loading Mash/2.3-GCC-12.3.0'
module load Mash/2.3-GCC-12.3.0

## copy sequences and split clusters
cp /PATH/all_plasmids.fasta ./output/all_plasmids.fasta
awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F="output/"ID".fasta"} {print >> F}' ./output/all_plasmids.fasta
rm ./output/all_plasmids.fasta

## loop through files
for file in /PATH/output/*.fasta;
do
        OUT=$(echo $file | cut -d / -f 4 | cut -d 'f' -f 1);
        mash screen -v 0.1 -i 0.99 -p ${SLURM_CPUS_PER_TASK} ${DB} ${file} > ./output/${OUT}tab
        # remove empty files, i.e. zero results
        find /PATH/output/ -type f -empty -delete
        rm ${file}
done

## Done
echo "Finished mash screen"
