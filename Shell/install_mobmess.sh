#!/bin/bash

#SBATCH --job-name=install_mobmess
#SBATCH --account=xx
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: install_mobmess.sh"

module purge

# load python
module load Python/3.12.3-GCCcore-13.3.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

folder=/PATH/mobmess

# create environment
echo "create python venv"
mkdir -p ${folder}
python -m venv ${folder}

# Activate the environment
echo "activate python venv"
source  ${folder}/bin/activate

# install dependencies (most are on Saga already)
pip install blosc
pip install pandas
pip install scipy
pip install scikit-learn

# Download and install PlasX
git clone https://github.com/michaelkyu/PlasX.git
pip install ./PlasX

# Download and install MobMess
git clone https://github.com/michaelkyu/MobMess.git
pip install ./MobMess

# deactivate venv
deactivate

# link programs
ln -s /PATH/delta-filter ${folder}/bin/delta-filter
ln -s /PATH/mummer ${folder}/bin/mummer
ln -s /PATH/nucmer ${folder}/bin/nucmer

echo "finished"
