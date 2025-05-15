#!/bin/bash

#SBATCH --job-name=geNomad
#SBATCH --account=xx
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: geNomad.sh"

module purge

# load python
module load Python/3.11.5-GCCcore-13.2.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

# Activate the environment
echo "activate python venv"
source  /PATH/bin/activate

# run geNomad
genomad end-to-end --threads ${SLURM_CPUS_PER_TASK} \
                   --cleanup /PATH/scaffolds.fasta \
                   PATH/output \
                   PATH/genomad_db

# deactivate venv
deactivate

echo "finished"
