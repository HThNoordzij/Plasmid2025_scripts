#!/bin/bash

#SBATCH --job-name=install_geNomad
#SBATCH --account=xx
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: install_geNomad.sh"

module purge

# load python
module load Python/3.11.5-GCCcore-13.2.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

folder=/PATH/genomad

# create environment
echo "create python venv"
mkdir -p ${folder}
python -m venv ${folder}

# Activate the environment
echo "activate python venv"
source  ${folder}/bin/activate

# install programs
pip install genomad

# download DB
genomad download-database .

# deactivate venv
deactivate

# link programs
ln -s /PATH/aragorn ${folder}/bin/aragorn
ln -s /PATH/mmseqs ${folder}/bin/mmseqs

echo "finished"
