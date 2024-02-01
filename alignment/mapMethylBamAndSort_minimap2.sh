#!/bin/bash
#SBATCH --job-name=mapMethylBamAndSort_minimap2.20231017
#SBATCH --partition=main
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=mapMethylBamAndSort_minimap2.20231017.%j.log
#SBATCH --time=18:00:00

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh; 
conda activate /private/home/mcechova/.conda/envs/methylation

in_cores=64 #number of processors to be used

unaligned_methyl_bam=$1
DIR="$(dirname "${unaligned_methyl_bam}")" #output in the same directory where the input file is
sample="$(basename ${unaligned_methyl_bam} .bam)"

ref_file=$2 #"/public/groups/migalab/mcechova/chm13v2.0.fa" #"GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref_name="$(basename -- $ref_file)"
ref_name="${ref_name%.*}"

in_args="-y -x map-ont --MD --eqx --cs -Y -L -a -k 17 -K 10g" #minimap parameters appropriate for nanopore reads

method="minimap2"

#script based on the wdl file by Melissa Meredith
#https://github.com/meredith705/ont_methylation/blob/32095600428d21bf53aef8a7ccc401b0f10a9145/tasks/minimap2.wdl

index_file="index.minimap2.${ref_name}.mmi"
if [ -f "$index_file" ]; then
    echo "$index_file exists."
else 
    echo "$index_file does not exist."
    #generate minimap index file
	minimap2 -k 17 -I 8G -d ${index_file} ${ref_file}
fi

#do the mapping with methylation tags
samtools fastq -T Mm,Ml ${unaligned_methyl_bam} | minimap2 -t ${in_cores} ${in_args} ${index_file} - | samtools view -@ ${in_cores} -bh - | samtools sort -@ ${in_cores} - > ${DIR}/${sample}.fastq.cpg.${method}.${ref_name}.bam
samtools index ${DIR}/${sample}.fastq.cpg.${method}.${ref_name}.bam

echo "Done."
