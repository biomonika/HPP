#!/bin/bash
#SBATCH --job-name=find_haplotype_for_flanks.sh.20240402
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --output=find_haplotype_for_flanks.sh.20240402.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=24
minIdentity=95

flanks=$1 #will be used to extract the flanks
flanks_name=$(basename -- "$flanks")
flanks_name="${flanks_name%.*}"

haplotype1=$2
haplotype1_name=$(basename -- "$haplotype1")
haplotype1_name="${haplotype1_name%.*}"

haplotype2=$3
haplotype2_name=$(basename -- "$haplotype2")
haplotype2_name="${haplotype2_name%.*}"

echo ${flanks_name} ${haplotype1_name} ${haplotype2_name}

chromosome=$(echo "$flanks" | cut -d'.' -f1)
order=$(echo "$flanks" | cut -d'.' -f8)

if [ -e "${flanks_name}.${haplotype1_name}" ]; then
    echo "Mashmap file with flanks mapped to haplotype 1 already exists."
else
	#map flanks to haplotype1
	mashmap --filter_mode one-to-one --threads ${threadCount} --perc_identity ${minIdentity} --noSplit --segLength 1000 -r ${haplotype1} -q ${flanks} -o ${flanks_name}.${haplotype1_name} > /dev/null 2>&1
fi

if [ -e "${flanks_name}.${haplotype2_name}" ]; then
    echo "Mashmap file with flanks mapped to haplotype 2 already exists."
else
	#map flanks to haplotype2
	mashmap --filter_mode one-to-one --threads ${threadCount} --perc_identity ${minIdentity} --noSplit --segLength 1000 -r ${haplotype2} -q ${flanks} -o ${flanks_name}.${haplotype2_name} > /dev/null 2>&1
fi

wait
# at this point, we mapped flanks to two haplotypes, and need to decide which one is a better match

identity_of_haplotype1_flank1=$(awk 'NR==1 {print $10}' "${flanks_name}.${haplotype1_name}")
identity_of_haplotype2_flank1=$(awk 'NR==1 {print $10}' "${flanks_name}.${haplotype2_name}")
identity_of_haplotype1_flank2=$(awk 'NR==2 {print $10}' "${flanks_name}.${haplotype1_name}")
identity_of_haplotype2_flank2=$(awk 'NR==2 {print $10}' "${flanks_name}.${haplotype2_name}")

echo "identity_of_haplotype1_flank1: " $identity_of_haplotype1_flank1
echo "identity_of_haplotype2_flank1: " $identity_of_haplotype2_flank1
echo "identity_of_haplotype1_flank2: " $identity_of_haplotype1_flank2
echo "identity_of_haplotype2_flank2: " $identity_of_haplotype2_flank2

# Compare decimal numbers using bc
if (( $(echo "$identity_of_haplotype1_flank1 > $identity_of_haplotype2_flank1" | bc -l) )) && \
   (( $(echo "$identity_of_haplotype1_flank2 > $identity_of_haplotype2_flank2" | bc -l) )); then
    echo "Haplotype 1 is a better match for both flanks."
elif (( $(echo "$identity_of_haplotype2_flank1 > $identity_of_haplotype1_flank1" | bc -l) )) && \
   (( $(echo "$identity_of_haplotype2_flank2 > $identity_of_haplotype1_flank2" | bc -l) )); then
    echo "Haplotype 2 is a better match for both flanks."
else
    echo "It is not possible to decide which haplotype is a better match."
fi

echo "Done."
echo "==========================="
date

