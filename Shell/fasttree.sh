#!/bin/bash

#SBATCH --job-name=FastTree
#SBATCH --account=xx
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

##Load modules
echo 'Loading FastTree/2.1.11-GCCcore-11.3.0'
module load FastTree/2.1.11-GCCcore-11.3.0

echo 'Usage: sbatch fasttree seq.fasta out'

##Run fasttree default settings
FastTree -nt -gtr -gamma $1 > $2

module purge
echo 'Done'
