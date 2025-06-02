#!/bin/bash

#SBATCH --job-name=bwa_samtools
#SBATCH --account=xx
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

## Variables
FASTA_FILE=$1
READ1=$2
READ2=$3
OUTPUT=$4

echo 'Usage: sbatch bwa_mapping.sh file.fasta read_R1.fasta read_R2.fasta output'

## Load modules
echo 'Loading BWA/0.7.17-GCCcore-11.3.0'
module load BWA/0.7.17-GCCcore-11.3.0

## Make index of FASTA_FILE
echo "indexing ${FASTA_FILE}"
bwa index -p IND_${OUTPUT} ${FASTA_FILE}

##echo "aligning reads to index with ${SLURM_CPUS_PER_TASK} threads"
bwa mem -t ${SLURM_CPUS_PER_TASK} IND_${OUTPUT} ${READ1} ${READ2} > output/${OUTPUT}.sam

## load modules
module purge
echo 'loading SAMtools/1.17-GCC-12.2.0'
module load SAMtools/1.17-GCC-12.2.0

### filtering and sorting reads
### with options -f 2 # reads paired in proper pair and -F 4 # removed unmapped reads
samtools view -@ ${SLURM_CPUS_PER_TASK} -S -b -F 4 -f 2 -e 'rlen>120' output/${OUTPUT}.sam | \
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o output/${OUTPUT}.sorted.bam -

##echo "index sorted bam file"
samtools index -@ ${SLURM_CPUS_PER_TASK} output/${OUTPUT}.sorted.bam

##echo "calculate coverage"
samtools coverage output/${OUTPUT}.sorted.bam -o output/${OUTPUT}_coverage.txt

## Done
echo "done"

## remove .sam file and IND files
rm output/${OUTPUT}.sam
rm IND_${OUTPUT}*
