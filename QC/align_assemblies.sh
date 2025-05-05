#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="align_assemblies.20250425"
#SBATCH --cpus-per-task=32
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --partition=long
#SBATCH --output=align_assemblies.20250425.%j.log

set -e
set -x

numCPU=32

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

# Input files and references
assembly1=$1
assembly2=$2

# Check if both files are provided
if [[ -z "$assembly1" || -z "$assembly2" ]]; then
  echo "Error: Two input files must be provided."
  exit 1
fi

# Check if files exist and are non-empty
if [[ ! -s "$assembly1" ]]; then
  echo "Error: File '$assembly1' does not exist or is empty."
  exit 1
fi

if [[ ! -s "$assembly2" ]]; then
  echo "Error: File '$assembly2' does not exist or is empty."
  exit 1
fi

echo "Both input files are valid. Proceeding..."

assembly1_name=$(basename "${assembly1%.*}")
assembly2_name=$(basename "${assembly2%.*}")

alignment_file="${assembly1_name}_${assembly2_name}.sam"
stat_file="${assembly1_name}_${assembly2_name}.stat.txt"
vcf_file="${assembly1_name}_${assembly2_name}.vcf"

minimap2 -t ${numCPU} --cs -x asm5 -c "$assembly1" "$assembly2" > "$alignment_file" 2> /dev/null
paftools.js stat "$alignment_file" >${stat_file}
paftools.js call -f "$assembly1" "$alignment_file" > "$vcf_file"

echo "Processing completed."
