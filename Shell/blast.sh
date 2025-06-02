#!/bin/bash

#SBATCH --job-name=blast_db
#SBATCH --account=xx
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

echo 'Usage: sbatch blast.sh'

##Load modules
echo 'Loading BLAST+/2.13.0-gompi-2022a'
module load BLAST+/2.13.0-gompi-2022a

## Run blastn
blastn  -query /PATH/all_plasmids.fasta \
        -task blastn \
        -db /PATH/db \
        -out all_plasmids_db.txt \
        -evalue 1 \
        -perc_identity 80 \
        -num_threads ${SLURM_CPUS_PER_TASK}\
        -outfmt "6 qseqid sseqid qlen slen qcovs qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"


## Message that you are done with the job
echo "Finished running"
