#!/bin/bash
#SBATCH --job-name=general.trf_analysis
#SBATCH --cpus-per-task=1
#SBATCH --array=1-184%46
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32g
#SBATCH --ntasks=1
#SBATCH --output=general.trf_analysis.20250126.%j.log
#SBATCH --time=18:00:00

set -e

pwd; hostname; date

# Specify the path to the config file
config=$1

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
assembly=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the assembly name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${assembly}"

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/microsatellite
#run TRF for mono and dinucleotides
trf ${assembly} 2 7 7 80 10 20 2 -d -h -l 23 || true #because trf returns the number of seqs as an exit code

echo "Done."
date
