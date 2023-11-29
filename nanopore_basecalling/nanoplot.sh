#!/bin/bash
#SBATCH --job-name=nanoplot.20231129
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=nanoplot.20231129.%j.log

pwd; hostname; date
set -e

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/nanoplot

fullfile=$1
file=$(basename -- "$fullfile")
file="${file%.*}"

numCPU=32

mkdir -p nanoplot.${file}
NanoPlot --threads ${numCPU} --fastq ${fullfile} --outdir nanoplot.${file} --prefix ${file}

echo "Done."
date

