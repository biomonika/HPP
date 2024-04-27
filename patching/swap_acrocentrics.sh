#!/bin/bash
#SBATCH --job-name=swap_acrocentrics.20240426
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=swap_acrocentrics.20240426.%j.log

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

#assembly_to_be_patched="PAN027.fully_phased.maternal.patched.fa"
#assembly_to_be_patched_mashmap="mashmap/PAN027.fully_phased.maternal.patched.mashmap.txt"
#assembly_to_be_patched_censat_track="censat/PAN027.fully_phased.acrocentric.MATERNAL.patched.cenSat.bed"
#assembly_reference="../maternal/duplex.maternal.scaffolds.fa"
#assembly_reference_mashmap="../maternal/mashmap/duplex.maternal.scaffolds.mashmap.txt"
#assembly_reference_censat_track="../maternal/censat/duplex.maternal.scaffolds.cenSat.bed"

assembly_to_be_patched="PAN027.fully_phased.paternal.patched.fa"
assembly_to_be_patched_mashmap="mashmap/PAN027.fully_phased.paternal.patched.mashmap.txt"
assembly_to_be_patched_censat_track="censat/PAN027.fully_phased.acrocentric.paternal.patched.cenSat.bed"
assembly_reference="../paternal/duplex.paternal.scaffolds.fa"
assembly_reference_mashmap="../paternal/mashmap/duplex.paternal.scaffolds.mashmap.txt"
assembly_reference_censat_track="../paternal/censat/duplex.paternal.scaffolds.cenSat.bed"

if [ ! -f "${assembly_to_be_patched}" ]; then
    echo "${assembly_to_be_patched} file not found. Exiting script."
    exit 1
fi
if [ ! -f "${assembly_to_be_patched_mashmap}" ]; then
    echo "${assembly_to_be_patched_mashmap} file not found. Exiting script."
    exit 1
fi
if [ ! -f "${assembly_to_be_patched_censat_track}" ]; then
    echo "${assembly_to_be_patched_censat_track} file not found. Exiting script."
    exit 1
fi
if [ ! -f "${assembly_reference}" ]; then
    echo "${assembly_reference} file not found. Exiting script."
    exit 1
fi
if [ ! -f "${assembly_reference_mashmap}" ]; then
    echo "${assembly_reference_mashmap} file not found. Exiting script."
    exit 1
fi
if [ ! -f "${assembly_reference_censat_track}" ]; then
    echo "${assembly_reference_censat_track} file not found. Exiting script."
    exit 1
fi


assembly_to_be_patched_name=$(basename -- "$assembly_to_be_patched")
assembly_to_be_patched_name="${assembly_to_be_patched_name%.*}"
assembly_reference_name=$(basename -- "$assembly_reference")
assembly_reference_name="${assembly_reference_name%.*}"

acrocentric_chromosome=("chr13" "chr14" "chr15" "chr21" "chr22")

swap_chromosome() {
	chromosome=$1

	#find where activeHOR coordinates start
	contig_name_assembly_to_be_patched=`grep ${chromosome} ${assembly_to_be_patched_mashmap} | cut -d' ' -f1`
	contig_length_assembly_to_be_patched=`grep ${chromosome} ${assembly_to_be_patched_mashmap} | cut -d' ' -f2`
	contig_name_assembly_reference=`grep ${chromosome} ${assembly_reference_mashmap} | cut -d' ' -f1`
	contig_length_assembly_reference=`grep ${chromosome} ${assembly_reference_mashmap} | cut -d' ' -f2`

	echo "contig_name_assembly_to_be_patched: $contig_name_assembly_to_be_patched"
	echo "contig_length_assembly_to_be_patched: $contig_length_assembly_to_be_patched"
	echo "contig_name_assembly_reference: $contig_name_assembly_reference"
	echo "contig_length_assembly_reference: $contig_length_assembly_reference"

	#we need to find activeHOR, and that's when we stop
	#this first part of the assembly will be ignored
	#even if we find multiple activeHOR hits, we are only interested in the first one of sufficient length
	cat ${assembly_to_be_patched_censat_track} | egrep ${contig_name_assembly_to_be_patched} | egrep "active_hor" | bedtools sort | bedtools merge | awk '$3 - $2 > 10000' >tmp.${chromosome}.activeHOR_start_assembly_to_be_patched.${contig_name_assembly_to_be_patched}.txt

	if [ ! -s "tmp.${chromosome}.activeHOR_start_assembly_to_be_patched.${contig_name_assembly_to_be_patched}.txt" ]; then
    	echo "Not sufficient activeHOR signal detected. Exiting script."
    	exit 1
	fi
	
	#write all hits into a temporary file, only use the first
	activeHOR_start_assembly_to_be_patched=`cat tmp.${chromosome}.activeHOR_start_assembly_to_be_patched.${contig_name_assembly_to_be_patched}.txt | head -n 1 | cut -f2`

	cat ${assembly_reference_censat_track} | egrep ${contig_name_assembly_reference} | egrep "active_hor" | bedtools sort | bedtools merge | awk '$3 - $2 > 10000' >tmp.${chromosome}.activeHOR_start_assembly_reference.${contig_name_assembly_reference}.txt

	if [ ! -s "tmp.${chromosome}.activeHOR_start_assembly_reference.${contig_name_assembly_reference}.txt" ]; then
    	echo "Not sufficient activeHOR signal detected. Exiting script."
    	exit 1
	fi

	#write all hits into a temporary file, only use the first
	activeHOR_start_assembly_reference=`cat tmp.${chromosome}.activeHOR_start_assembly_reference.${contig_name_assembly_reference}.txt | head -n 1 | cut -f2`

	#remove unnecessary temporary files
	rm -f tmp.${chromosome}.activeHOR_start_assembly_to_be_patched.${contig_name_assembly_to_be_patched}.txt
	rm -f tmp.${chromosome}.activeHOR_start_assembly_reference.${contig_name_assembly_reference}.txt

	#only keep part of the assembly AFTER the first activeHOR
	region=${contig_name_assembly_to_be_patched}:${activeHOR_start_assembly_to_be_patched}-${contig_length}
	echo ${region}
	samtools faidx ${assembly_to_be_patched} ${region} >${chromosome}.${assembly_to_be_patched_name}.toKeep.fa


	#only keep part of the assembly BEFORE the first activeHOR
	region=${contig_name_assembly_reference}:0-${activeHOR_start_assembly_reference}
	echo ${region}
	samtools faidx ${assembly_reference} ${region} >${chromosome}.${assembly_reference_name}.toKeep.fa

	#COMBINE BOTH ASSEMBLIES
	#add header
	echo ">${chromosome}.${assembly_to_be_patched_name}.acroswap" >tmp.${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa
	#add sequence
	cat ${chromosome}.${assembly_reference_name}.toKeep.fa ${chromosome}.${assembly_to_be_patched_name}.toKeep.fa | seqtk seq | egrep -v "^>" | tr -d '\n' >>tmp.${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa

	#reformat to 60 characters per line in fasta file
	seqtk seq -l 60 tmp.${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa >${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa
	rm tmp.${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa

	#remove unnecessary files
	ls ${chromosome}.${assembly_to_be_patched_name}.toKeep.fa
	ls ${chromosome}.${assembly_reference_name}.toKeep.fa
	rm ${chromosome}.${assembly_to_be_patched_name}.toKeep.fa ${chromosome}.${assembly_reference_name}.toKeep.fa

	echo "After acroswap, the resulting patched file is: " ${chromosome}.PATCHED.${assembly_to_be_patched_name}.acroswap.fa.acroswap.fa
}

for chromosome in "${acrocentric_chromosome[@]}"; do
    echo "$chromosome"
    swap_chromosome $chromosome
done

echo "Done."
date
