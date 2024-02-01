#!/bin/bash

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh; 
conda activate /private/home/mcechova/.conda/envs/methylation

in_cores=200 #number of processors to be used

unaligned_methyl_bam=$1
DIR="$(dirname "${unaligned_methyl_bam}")" #output in the same directory where the input file is
sample="$(basename ${unaligned_methyl_bam} .bam)"

ref_file=$2 #"/public/groups/migalab/mcechova/chm13v2.0.fa" #"GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref_name="$(basename -- $ref_file)"
ref_name="${ref_name%.*}"

method="winnowmap"

#check if appropriate meryl database exists yet

meryl=repetitive_${ref_name}_k15.txt
if [ -f "$FILE" ]; then
    echo "${meryl} database already exists."
else 
    echo "${meryl} database needs to be created."
    meryl count k=15 output merylDB_${ref_name}_k15 ${ref_file}
	meryl print greater-than distinct=0.9998 merylDB_${ref_name}_k15 > ${meryl}
fi


#script based on the wdl file by Melissa Meredith
#https://github.com/meredith705/ont_methylation/blob/32095600428d21bf53aef8a7ccc401b0f10a9145/tasks/minimap2.wdl

samtools fastq -T MM,ML ${unaligned_methyl_bam} | winnowmap --cs -t ${in_cores} -I 8G -W ${meryl} -ax map-ont -y ${ref_file} - | samtools view -@ ${in_cores} -bh - | samtools sort -@ ${in_cores} - > ${DIR}/${sample}.fastq.cpg.${method}.${ref_name}.bam
samtools index ${DIR}/${sample}.fastq.cpg.${method}.${ref_name}.bam

echo "Done."

