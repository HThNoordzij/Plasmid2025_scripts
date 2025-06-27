#!/bin/bash

#SBATCH --job-name=kalign
#SBATCH --account=xx
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

##Load modules
echo 'Loading Kalign/3.3.5-GCCcore-11.3.0'
module load Kalign/3.3.5-GCCcore-11.3.0

echo 'Usage: sbatch kalign.sh seq.fasta out'

##Set variables

##Run kalign with default settings
kalign -i $1 \
       -o $2.fasta \
       --nthreads ${SLURM_CPUS_PER_TASK}

module purge
echo 'Done'
