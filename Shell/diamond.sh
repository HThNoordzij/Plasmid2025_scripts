#!/bin/bash

#SBATCH --job-name=diamond
#SBATCH --account=xx
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

## Variables
FASTA_FILE=$1
DATABASE_NAME=$2
OUTPUT=$3
echo 'Usage: sbatch diamond.sh file.fasta database output_name'

##Load modules
echo 'Loading DIAMOND/2.1.0-GCC-11.3.0'
module load DIAMOND/2.1.0-GCC-11.3.0

## Run diamond
echo  "Running diamond on ${FASTA_FILE} with database ${DATABASE_NAME}. Output file: ${OUTPUT}_diamond_out6.txt"
diamond blastp -d ${DATABASE_NAME} \
        -q ${FASTA_FILE} \
        -o ${OUTPUT}_diamond_outfmt6.tsv \
        --threads ${SLURM_CPUS_PER_TASK}\
        --outfmt 6 qseqid sseqid qlen slen qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

## Message that you are done with the job
echo "Finished running"
