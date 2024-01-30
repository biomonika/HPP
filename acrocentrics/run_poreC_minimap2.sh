#!/bin/bash
#SBATCH --job-name=run_poreC_minimap2.sh.20240121
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=run_poreC_minimap2.sh.20240121.%j.log

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

reference=$1 #${individual}.diploid.complete.fa
porec_file=$2 #${individual}.poreC.fastq.gz
filename=$(basename "$reference" .fa)."poreC"

threadCount=32
minimap2_args="-a -x map-ont -k 17 -K 10g -I 8g"
min_mapq=10

#parameters from https://github.com/meredith705/gfase_wdl/blob/main/tasks/proximityAlign.wdl

# align with minimap2
minimap2 ${minimap2_args} \
         -t ${threadCount} \
         ${reference} \
         ${porec_file} | samtools view -bh -@ ${threadCount} -q ${min_mapq} -o ${filename}.bam -O BAM -

samtools view ${filename}.bam | cut -f1,3,4 | sort -k1,1 -k3,3n >${filename}.txt
input_minimap2=${filename}.txt 

#create all pairwise combination for poreC reads
python combinations.py ${input_minimap2}

#select relevant columns
cat output_${input_minimap2} | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$4}' >pairwise.minimap2.${filename}.complete.txt

#summarize all interactions
cat pairwise.minimap2.${filename}.complete.txt | cut -f1,4 | sort | uniq -c | sort -rgk1 | tr -s ' ' | sed s'/^ //g' | sed s'/\t/ /g' >summary.pairwise.minimap2.${filename}.complete.txt

rm pairwise.minimap2.${filename}.complete.txt


echo "Done."
date

