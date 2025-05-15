#!/bin/bash

#SBATCH --job-name=mobmess
#SBATCH --account=xx
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: mobmess.sh"

workfolder=/PATH/mobmess

module purge

# load python
module load Python/3.12.3-GCCcore-13.3.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

# Activate the environment
echo "activate python venv"
source  /PATH/bin/activate

mkdir ${workfolder}

# run mobmess
mobmess systems \
    --sequences /PATH/scaffolds_geNomad.fasta \
    --complete /PATH/circular_scaffolds_geNomad.txt \
    --output /PATH/output \
    --threads  ${SLURM_CPUS_PER_TASK} \
    --tmp ${workfolder}

# deactivate venv
deactivate

echo "finished"
