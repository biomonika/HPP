#!/bin/bash

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /public/home/mcechova/conda/hifiasm


in_cores=200 #number of processors to be used

mapped_methyl_bam=$1 #provide the diploid alignment
INPUT_DIR="$(dirname "${mapped_methyl_bam}")" #the same directory where the input file is
OUTPUT_DIR="/data/mcechova/CDR/secphase" #where the output files should be stored ###TO BE MODIFIED###
cd ${OUTPUT_DIR}
BAM_PREFIX="$(basename ${mapped_methyl_bam} .bam)"

ref_file=$2 #the reference used for the diploid alignment
ref_name="$(basename -- $ref_file)"
FASTA_PREFIX="${ref_name%.*}"

## Sort by read name
samtools sort -n -@${in_cores} ${mapped_methyl_bam} > ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam

## Phase reads for ONT
docker run \
    -v ${INPUT_DIR}:${INPUT_DIR} \
    mobinasri/secphase:v0.1 \
    secphase -q -c -t20 -d 1e-2 -e 0.1 -b20 -m10 -s20 \
    -i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
    -f ${INPUT_DIR}/${FASTA_PREFIX}.fa > ${INPUT_DIR}/${BAM_PREFIX}.secphase.log

echo "Phasing finished."

#correct
docker run \
    -v ${INPUT_DIR}:${INPUT_DIR} \
    -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
    mobinasri/secphase:v0.1 \
    correct_bam \
    -i ${INPUT_DIR}/${BAM_PREFIX}.bam \
    -P ${INPUT_DIR}/${BAM_PREFIX}.secphase.log \
    -o ${OUTPUT_DIR}/${BAM_PREFIX}.corrected.bam \
    --primaryOnly

echo "Correction finished."
echo "Done."
