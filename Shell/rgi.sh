#!/usr/bin/env bash

#SBATCH --job-name=rgi_plasmids
#SBATCH --account=xx
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

module purge

# Load db
singularity run --bind /PATH/ /PATH/rgi/rgi.sif rgi load --card_json /PATH/card_3.3.0/card.json --local # version 3.3.0

# Run rgi
singularity run --bind /PATH/ /PATH/rgi/rgi.sif rgi main \
        --input_sequence /PATH/scaffolds/IDx_plasmids.fasta \
        --output_file /PATH/rgi/output/idx \
        --local \
        --clean \
        --keep \
        --low_quality \
        --num_threads ${SLURM_CPUS_PER_TASK} \
        --split_prodigal_jobs
