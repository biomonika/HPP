#!/bin/bash
#SBATCH --job-name=QC_trim.20231127
#SBATCH --partition=main
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=QC_trim.20231127.%j.log
#SBATCH --time=3:00:00

pwd; hostname; date
set -e

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/ONT

fullfile=$1
file=$(basename -- "$fullfile")
file="${file%.*}"

numCPU=32

#TRIM AND FILTER THE DUPLEX READS
seqkit subseq --threads ${numCPU} -r 100:-10 ${fullfile}| chopper --threads ${numCPU} -q 30 -l 15000 | pigz -p ${numCPU} > ${file}.Q30.15kb.fastq.gz
fastqc --threads ${numCPU} ${file}.Q30.15kb.fastq.gz

seqkit subseq --threads ${numCPU} -r 100:-10 ${fullfile}| chopper --threads ${numCPU} -q 25 -l 15000 | pigz -p ${numCPU} > ${file}.Q25.15kb.fastq.gz
fastqc --threads ${numCPU} ${file}.Q25.15kb.fastq.gz

#GENERATE QC PLOTS
conda activate /private/home/mcechova/conda/nanoplot
mkdir -p nanoplot.${file}

NanoPlot --threads ${numCPU} --fastq ${file}.Q30.15kb.fastq.gz --outdir nanoplot.${file} --prefix Q30.15kb
NanoPlot --threads ${numCPU} --fastq ${file}.Q25.15kb.fastq.gz --outdir nanoplot.${file} --prefix Q25.15kb

echo "Done."
date

