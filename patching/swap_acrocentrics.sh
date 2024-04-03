#!/bin/bash
#SBATCH --job-name=swap_acrocentrics.20240327
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=swap_acrocentrics.20240327.%j.log

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

#I need the following files:
#1. assembly that needs patching
#2. associated mashmap files of the assembly that needs patching
#3. associated censat track of the assembly that needs patching
#4. assembly that will be used for patching
#5. associated mashmap files of the assembly that will be used for patching
#6. associated censat track of the assembly that will be used for patching
#7. the name of the chromosome that we want to patch

assembly_to_be_patched=$1
assembly_to_be_patched_mashmap=$2
assembly_to_be_patched_censat_track=$3
assembly_reference=$4
assembly_reference_mashmap=$5
assembly_reference_censat_track=$6
chromosome=$7 #which chromosome are we trying to patch

assembly_to_be_patched="full.maternal.contigs.fa"
assembly_to_be_patched_mashmap="mashmap/full.maternal.contigs.mashmap.txt"
assembly_to_be_patched_censat_track="full.maternal.contigs.cenSat.bed"
assembly_reference="duplex.maternal.contigs.fa"
assembly_reference_mashmap="mashmap/duplex.maternal.contigs.mashmap.txt"
assembly_reference_censat_track="duplex.maternal.contigs.cenSat.bed"

assembly_to_be_patched_name=$(basename -- "$assembly_to_be_patched")
assembly_to_be_patched_name="${assembly_to_be_patched_name%.*}"
assembly_reference_name=$(basename -- "$assembly_reference")
assembly_reference_name="${assembly_reference_name%.*}"

acrocentric_chromosome=("chr13" "chr14" "chr15" "chr21" "chr22")

swap_chromosome() {
	chromosome=$1

	#find where rDNA coordinates start
	contig_name_assembly_to_be_patched=`grep ${chromosome} ${assembly_to_be_patched_mashmap} | cut -d' ' -f1`
	contig_length_assembly_to_be_patched=`grep ${chromosome} ${assembly_to_be_patched_mashmap} | cut -d' ' -f2`
	contig_name_assembly_reference=`grep ${chromosome} ${assembly_reference_mashmap} | cut -d' ' -f1`
	contig_length_assembly_reference=`grep ${chromosome} ${assembly_reference_mashmap} | cut -d' ' -f2`

	#we need to find rDNA of length at least 1kb, and that's when we stop
	#this first part of the assembly will be ignored
	#even if we find multiple rDNA hits, we are only interested in the first one of sufficient length
	rDNA_start_assembly_to_be_patched=`cat ${assembly_to_be_patched_censat_track} | egrep ${contig_name_assembly_to_be_patched} | egrep "rDNA" | bedtools sort | bedtools merge | awk '$3 - $2 > 5000' | head -n 1 | cut -f2 `
	rDNA_start_assembly_reference=`cat ${assembly_reference_censat_track} | egrep ${contig_name_assembly_reference} | egrep "rDNA" | bedtools sort | bedtools merge | awk '$3 - $2 > 5000' | head -n 1 | cut -f2 `

	#only keep part of the assembly AFTER the first rDNA
	region=${contig_name_assembly_to_be_patched}:${rDNA_start_assembly_to_be_patched}-${contig_length}
	echo ${region}
	samtools faidx ${assembly_to_be_patched} ${region} >${assembly_to_be_patched_name}.toKeep.fa


	#only keep part of the assembly BEFORE the first rDNA
	region=${contig_name_assembly_reference}:0-${rDNA_start_assembly_reference}
	echo ${region}
	samtools faidx ${assembly_reference} ${region} >${assembly_reference_name}.toKeep.fa

	#COMBINE BOTH ASSEMBLIES
	#add header
	echo ">${chromosome}.${assembly_to_be_patched_name}.acroswap" >PATCHED.${chromosome}.${assembly_to_be_patched_name}.acroswap.fa
	#add sequence
	cat ${assembly_reference_name}.toKeep.fa ${assembly_to_be_patched_name}.toKeep.fa | seqtk seq | egrep -v "^>" | tr -d '\n' >>PATCHED.${chromosome}.${assembly_to_be_patched_name}.acroswap.fa

	#remove unnecessary files
	ls ${assembly_to_be_patched_name}.toKeep.fa
	ls ${assembly_reference_name}.toKeep.fa
	rm ${assembly_to_be_patched_name}.toKeep.fa ${assembly_reference_name}.toKeep.fa

	echo "After acroswap, the resulting patched file is: " PATCHED.${chromosome}.${assembly_to_be_patched_name}.acroswap.fa.acroswap.fa
}

for chromosome in "${acrocentric_chromosome[@]}"; do
    echo "$chromosome"
    swap_chromosome $chromosome
done

echo "Done."
date

