#!/bin/bash
#SBATCH --job-name=general.microsatellite_analysis
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8%8
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32g
#SBATCH --ntasks=8
#SBATCH --output=general.microsatellite_analysis.20250126.%j.log
#SBATCH --time=24:00:00

set -e

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/QC

# Specify the path to the config file
config=$1

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
assembly=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${assembly}"

#run seqrequester
/private/home/mcechova/seqrequester/build/bin/seqrequester microsatellite -prefix ${assembly_name}.microsatellite -window 128 -ga -gc -at ${assembly}
#remove unnecessary raw files
rm -f *microsatellite.AT.bed
rm -f *microsatellite.GA.bed
rm -f *microsatellite.GC.bed
rm -f *microsatellite.TC.bed

#run homopolymer annotation by Nancy
python make_mononuc_bedfile.py --fasta ${assembly} --minsize 10 >${assembly_name}.mono.bed

conda activate /private/home/mcechova/conda/microsatellite
#run TRF for mono and dinucleotides
trf ${assembly} 2 7 7 80 10 50 2 -d

echo "Done."
date
