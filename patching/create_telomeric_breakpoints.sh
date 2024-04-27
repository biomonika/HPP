#!/bin/bash
#SBATCH --job-name=create_telomeric_breakpoints.20240427
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=create_telomeric_breakpoints.20240427.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/ONT

assembly=$1
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"
#assembly="full.maternal.contigs.fa"
flank_size=1000000 #set the default flank length to 1Mbp

#store all chromosome names in a bash array
chromosomes=()
while read -r chromosome; do
    # Add the chromosome name to the array
    chromosomes+=("$chromosome")
done < <(bioawk -c fastx '{print $name}' "$assembly")

seqtk telo ${assembly} >${assembly_name}.telomeres.txt

#create the files to which we want to store the information about telomeres at the beginnings and ends of chromosomes
touch ${assembly_name}.telomeres.start.txt
touch ${assembly_name}.telomeres.end.txt

while IFS= read -r line; do
    echo $line

    start=""
    end=""
    length=""

    # Split the line into columns
    IFS=$'\t' read -r chromosome start end length <<< "$line"
    echo "start: $start" 
    echo "end: $end" 
    echo "length: $length" 
    
    # Telomere is present at the beginning
    if [ "$start" -eq 0 ]; then
        echo "Start coordinate is 0 for chromosome: $chromosome"
        echo $line >${assembly_name}.telomeres.start.txt
    fi

    # Telomere is present at the end
    if [ "$end" -eq "$length" ]; then
        echo "End coordinate is equal to length for chromosome: $chromosome"
        echo $line >${assembly_name}.telomeres.end.txt
    fi
done < "${assembly_name}.telomeres.txt"


#now for each chromosome, we check if telomere is present at both start and end
for chromosome in "${chromosomes[@]}"; do
    echo "Chromosome: $chromosome"

    if grep -q "$chromosome" "${assembly_name}.telomeres.start.txt"; then
        echo "Telomere found in ${assembly_name}.telomeres.start.txt"
    else
        echo "Telomere NOT found in ${assembly_name}.telomeres.start.txt"
        echo "Bed file will be created."
        echo -e "$chromosome\t0\t$flank_size" > "${chromosome}.${assembly_name}.telomeres.start.bed"
    fi

    if grep -q "$chromosome" "${assembly_name}.telomeres.end.txt"; then
        echo "Telomere found in ${assembly_name}.telomeres.end.txt"
    else
        echo "Telomere NOT found in ${assembly_name}.telomeres.end.txt"
        echo "Bed file will be created."
        samtools faidx ${assembly} $chromosome >tmp.${chromosome}.fa
        chromosome_length=$(bioawk -c fastx '{ print length($seq) }' < "tmp.${chromosome}.fa")
        rm -r "tmp.${chromosome}.fa"
        flank_threshold=$((chromosome_length - flank_size))
        echo -e "$chromosome\t$flank_threshold\t$chromosome_length" > "${chromosome}.${chromosome_length}.telomeres.start.bed"
        bedtools getfasta -fi ${assembly} -bed "${chromosome}.${chromosome_length}.telomeres.start.bed" >"${chromosome}.${chromosome_length}.telomeres.start.bed.fa"

    fi
done

#remove unnecessary files
rm -r ${assembly_name}.telomeres.start.txt
rm -r ${assembly_name}.telomeres.end.txt

echo "Done."
date

