#!/bin/bash

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /public/home/mcechova/conda/alignment

myThreads=50

cd "/data/mcechova/guppy"

extracted_raw_folder=$1 #provide a path to the extracted tar file
tar_file_name="$(basename -- $extracted_raw_folder)"
basecalled_folder="/data/mcechova/guppy/basecalled_5hmc_5mc_${tar_file_name}"

#check if basecalled files already exist

if [ -f ${basecalled_folder}/basecalled_${tar_file_name}.pass.bam} ] 
then
	echo "Basecalled and merged files already exist. Skipping basecalling and merging." 
else
	echo "Basecalled and merged files DO NOT exist."
	find ${extracted_raw_folder} -type f -name '*.fast5' | wc -l

	guppy_basecaller \
		-i ${extracted_raw_folder} \
		-s ${basecalled_folder} \
		-c /opt/ont/guppy/data/dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_sup_prom.cfg \
		--bam_out \
		--index \
		--compress_fastq \
		-x cuda:4,5,6,7 \
		-r \
		--read_batch_size 25000 \
		--records_per_fastq 25000

	#remove the extracted raw data to save space once the basecalling gas finished
	rm -r ${extracted_raw_folder}

	#merge bam file
	samtools merge --threads ${myThreads} ${basecalled_folder}/basecalled_${tar_file_name}.pass.bam ${basecalled_folder}/pass/*.bam
	samtools merge --threads ${myThreads} ${basecalled_folder}/basecalled_${tar_file_name}.fail.bam ${basecalled_folder}/fail/*.bam

	#merge fastq.gz
	cat ${basecalled_folder}/pass/*.bam >${basecalled_folder}/basecalled_${tar_file_name}.pass.fastq.gz
	cat ${basecalled_folder}/fail/*.bam >${basecalled_folder}/basecalled_${tar_file_name}.fail.fastq.gz

	#remove unneccessary folder
	rm -r ${basecalled_folder}/pass
	rm -r ${basecalled_folder}/fail

	samtools index ${basecalled_folder}/basecalled_${tar_file_name}.pass.bam
	samtools index ${basecalled_folder}/basecalled_${tar_file_name}.fail.bam
fi

samtools flagstat ${basecalled_folder}/basecalled_${tar_file_name}.pass.bam
samtools flagstat ${basecalled_folder}/basecalled_${tar_file_name}.fail.bam

echo "Done."




