#!/bin/bash
#SBATCH --job-name=run_HiC_bwa_mem.sh.20240118
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=24
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=run_HiC_bwa_mem.sh.20240118.%j.log


set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=24

#if using this code, the input files are expected to be derived from the name of an individual, for example:
#PAN027.R1.fastq.gz
#PAN027.R2.fastq.gz
#PAN027.diploid.complete.fa

individual=$1
echo ${individual}

#USE ILLUMINA READS AND THE REFERENCE AS THE INPUTS
hic_file="${individual}.R1.fastq.gz ${individual}.R2.fastq.gz" #Hi-C reads
reference=${individual}.diploid.complete.fa #reference

min_mapq=1

filename=$(basename "$reference" .fa)."HiC"

bwa index ${reference}

# align with bwa mem since using Illumina reads
bwa mem \
         -t ${threadCount} \
         ${reference} \
         ${hic_file} | samtools view -bh -@ ${threadCount} -q ${min_mapq} -o ${filename}.HiC.bam -O BAM -

samtools view ${filename}.HiC.bam | cut -f1,3,4 | sort -k1,1 -k3,3n >${filename}.txt
input_bwamem=${filename}.txt 
python combinations.py ${input_bwamem}

#select relevant columns
cat output_${input_bwamem} | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$4}' >pairwise.bwamem.${filename}.complete.txt

#summarize all interactions
cat pairwise.bwamem.${filename}.complete.txt | cut -f1,4 | sort | uniq -c | sort -rgk1 | tr -s ' ' | sed s'/^ //g' | sed s'/\t/ /g' >summary.pairwise.bwamem.${filename}.complete.txt

rm pairwise.bwamem.${filename}.complete.txt

echo "Done."
date

