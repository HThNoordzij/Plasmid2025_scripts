#!/bin/bash

#SBATCH --job-name=macsyfinder_install
#SBATCH --account=xx
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: install_macsyfinder.sh"

module purge

# load python
module load Python/3.12.3-GCCcore-13.3.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

# create environment
echo "create python venv"
python -m venv /PATH/venv/macsyfinder/


# Activate the environment
echo "activate python venv"
source /PATH/venv/macsyfinder/bin/activate

# install macsyfinder
pip install macsyfinder
pip install hmmer

# install models
macsydata install CONJScan
macsydata install TFFscan
macsydata install TXSScan
macsydata install CasFinder

# deactivate venv
deactivate

echo "finished"
