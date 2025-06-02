#!/bin/bash

#SBATCH --job-name=kraken2
#SBATCH --account=xx
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

## Variables
DB=/PATH/
FASTA_FILE=$1
OUT_FILE=$2

echo 'Usage: sbatch kraken2.sh fasta.file outname'

## Load modules
echo 'Kraken2/2.1.2-gompi-2022a'
module load Kraken2/2.1.2-gompi-2022a

## run kraken
kraken2 --db ${DB} \
        --threads ${SLURM_CPUS_PER_TASK} \
        --report ${OUT_FILE}.txt \
        --use-names \
        --classified-out ${OUT_FILE}.fasta \
        --log ${OUT_FILE}_log.txt \
        ${FASTA_FILE}
