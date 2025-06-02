#!/bin/bash

#SBATCH --job-name=macsyfinder_CONJscan
#SBATCH --account=xx
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variable as an error

echo "Usage: macsyfinder.sh"

module purge

# load python
module load Python/3.12.3-GCCcore-13.3.0

# Set the ${PS1} (needed in the source of the virtual environment for some Python versions)
export PS1=\$

# Activate the environment
echo "activate python venv"
source /PATH/venv/macsyfinder/bin/activate

# loop over ORF files per plasmid and run macsyfinder
for file in /PATH/*.fasta;
do
        BASE=$(echo $file | cut -d / -f 3);  # -f 3 is depending on the PATH length
        OUT=$(echo $BASE | cut -d '_' -f 1);

        macsyfinder -m CONJScan/Plasmids all \
                --sequence-db $file \
                -w ${SLURM_CPUS_PER_TASK} \
                --db-type ordered_replicon \
                --replicon-topology circular \
                -o /PATH/$OUT
done

# deactivate venv
deactivate

echo "finished"
