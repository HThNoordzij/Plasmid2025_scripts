#!/bin/bash

#SBATCH --job-name=metaplasmidSPAdes
#SBATCH --account=xx
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=10G

##Usage: sbatch metaplasmidSpades.sh R1.fasta R2.fasta output/

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

##Load modules
echo 'Loading SPAdes'
module load SPAdes/3.15.3-GCC-10.3.0

##Set variables
READ1=$1
READ2=$2
OUTPUT=$3
let MEM_R=250

echo 'Usage: sbatch metaplasmidSPAdes.sh R1.fasta R2.fasta output'

##Run metaplasmidspades
metaplasmidspades.py \
	-1 ${READ1} \
	-2 ${READ2} \
	-o ${OUTPUT} \
	--threads ${SLURM_CPUS_PER_TASK} \
	--memory ${MEM_R} 

echo 'Done'
