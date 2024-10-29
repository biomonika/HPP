#!/bin/bash
#SBATCH --job-name=bedGraphToBigWig
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=bedGraphToBigWig.20241029.%j.log
#SBATCH --time=23:00:00

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/groups/migalab/jmmenend/.conda_envs/ucsc

# Check if the script receives exactly 3 arguments
if [ "$#" -ne 3 ]; then
  echo "Error: Three arguments are required."
  echo "Usage: $0 <bedgraph> <chrom_sizes> <output_file>"
  exit 1
fi

# Assign input parameters to variables
bedgraph=$1
chrom_sizes=$2
output_file=$3

# Check if the bedgraph file has at least 4 columns
column_count=$(awk '{print NF; exit}' "$bedgraph")
if [ "$column_count" -lt 4 ]; then
  echo "Error: The bedgraph file must contain at least 4 columns."
  exit 1
fi

if [ "$column_count" -eq 4 ]; then

  bedGraphToBigWig "$bedgraph" "$chrom_sizes" "$output_file"

elif [ "$column_count" -gt 4 ]; then

  # Subset to only the first 4 columns if there are more
  awk '{print $1, $2, $3, $4}' "$bedgraph" > "${bedgraph}.tmp"
  bedGraphToBigWig "${bedgraph}.tmp" "$chrom_sizes" "$output_file"
  rm "${bedgraph}.tmp"

fi

echo "Done."
