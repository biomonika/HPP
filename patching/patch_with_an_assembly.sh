#!/bin/bash
#SBATCH --job-name=patch_with_an_assembly.20240325
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --output=patch_with_an_assembly.20240325.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=24
minIdentity=95
reference="chm13v2.0.fa.gz"
bed_file=$1
assembly=$2
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"
patch_reference=$3
patch_reference_name=$(basename -- "$patch_reference")
patch_reference_name="${patch_reference_name%.*}"

echo ${bed_file} ${patch_reference} ${haplotype}

chromosome=$(echo "$bed_file" | cut -d'.' -f1)
order=$(echo "$bed_file" | cut -d'.' -f8)

#extract flanks and find out where they belong
bedtools getfasta -fi ${assembly} -bed ${bed_file} -name >flanks.${bed_file}.${patch_reference_name}

#BREAKPOINTS TO AN ASSEMBLY AVAILABLE FOR PATCHING

if [ -e "${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file exists and won't be re-written."
else
    echo "Wfmash does not exist. Creating now."
    wfmash --threads ${threadCount} --segment-length=1000 --map-pct-id=${minIdentity} --no-split ${patch_reference} ${flank_file} >tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    cat tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt | sed s'/\t/ /g' | cut -d' ' -f1-10 | sort -k1,1n >${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    #remove temporary wfmash file
    rm tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
fi

wait

file_names=("${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt")

max_avg=0
max_file=""

# Loop through each file
for file_name in "${file_names[@]}"; do
    echo "Analyzing mashmap file with flanks mapped to the reference used for patching. "
    num_lines=$(wc -l < "$file_name")

    # Check if the number of lines is equal to two
    if [ "$num_lines" -eq 2 ]; then
        # Check if the 6th column is identical within the file
        if [ "$(awk '{print $6}' "$file_name" | uniq | wc -l)" -eq 1 ]; then
            
            first_coordinate_of_gap=$(awk 'NR==1 {print $9}' "${file_name}")
            second_coordinate_of_gap=$(awk 'NR==2 {print $8}' "${file_name}")
            test_gap=$((second_coordinate_of_gap - first_coordinate_of_gap))
        
            echo "--------"
            echo "$file_name"
            echo "test_gap $test_gap"
            echo "--------"

            # Calculate the average of numbers in the 10th column
            avg=$(awk '{sum += $NF} END {print sum/NR}' "$file_name")

            # Compare averages
            if (( $(echo "$avg > $max_avg" | bc -l) )); then
                max_avg=$avg
                max_file=$file_name
            fi
        else
            echo "Skipping $file_name as the 6th column is not identical (different contig names)."
            echo "The breakpoint is not resolved in the assembly used for patching."
        fi
    else
        echo "Skipping $file_name as the file does not have two rows. "
        echo "One of the flanks is not mapped successfully."
    fi
done

echo ""

if [ -n "$max_file" ]; then
    echo "The file with the highest average in the 10th column and identical 6th column is: $max_file"
else
    echo "No valid file found with identical 6th column."
    rm -f flanks.${bed_file}.${patch_reference_name}
    exit -1
fi

echo ""

#Get coordinates of the GAP/BREAKPOINT REGION
gap_start=""
gap_end=""
contig_name=""
extract_variables() {
    file_name=$1
    # Extract 9th column from the first row
    gap_start=$(awk 'NR==1 {print $9}' "$file_name")
    # Extract 10th column from the second row
    gap_end=$(awk 'NR==2 {print $8}' "$file_name")
    # Extract contig name
    contig_name=$(awk 'NR==1 {print $6}' "$file_name")

    # Print the extracted variables
    echo "For file $file_name:"
    echo "GAP/BREAKPOINT contig name (6th column): $contig_name"
    echo "GAP/BREAKPOINT start (9th column from the first row): $gap_start"
    echo "GAP/BREAKPOINT end (8th column from the second row): $gap_end"
}

# Extract coordinates from the reference used for patching
extract_variables "$max_file"

gap_size=$((gap_end - gap_start))
echo "Gap size is $gap_size"

region=""
# Check if the gap size is negative
if [ "$gap_size" -lt 0 ]; then
    echo "The gap_size is negative. We will be extracting reverse complement of the sequence."
    region=${contig_name}:${gap_end}-${gap_start}
    samtools faidx ${patch_reference} ${region} >to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
    seqtk seq -r to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa >${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
    rm to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
else
    echo "The gap_size is positive."
    region=${contig_name}:${gap_start}-${gap_end}
    samtools faidx ${patch_reference} ${region} >${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
fi

# Extract sequence names of the contigs we will be merging
echo "Extract sequence names of the contigs we will be merging"
chr_mashmap=${bed_file}
sequence_names=$(awk '{print $1}' ${chr_mashmap})

first_contig=$(awk 'NR==1 {print $1}' "${chr_mashmap}")
second_contig=$(awk 'NR==2 {print $1}' "${chr_mashmap}")

while IFS= read -r sequence_name; do
    # Use samtools faidx to extract the sequence
    # We are extracting the original sequences from the assembly
    echo $sequence_name
    samtools faidx ${assembly} ${sequence_name} > ${sequence_name}.fasta
done <<< "$sequence_names"

# Use head to extract the first line of each file
header1=$(head -n 1 "${first_contig}.fasta")
header2=$(head -n 1 "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa")
header3=$(head -n 1 "${second_contig}.fasta")

#MERGE SEPARATELY HEADER AND THE BODY/SEQUENCE
echo ">"${chromosome}.$(echo "original_$header1" | tr -d '>')"+"$(echo "${patch_reference_name}.${gap_size}.${header2}" | tr -d '>')"+"$(echo "original_$header3" | tr -d '>') >${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "${first_contig}.fasta" | grep -v ">" >>${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa" | grep -v ">" >>${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "${second_contig}.fasta" | grep -v ">" >>${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp

#reformat the final patched fasta chromosome
seqtk seq -L 0 ${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp >${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta
rm -f ${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta.tmp
echo "Merged concatenated fasta for ${chromosome} saved to" ${chromosome}.${order}.PATCHED.${assembly_name}.with.${patch_reference_name}.fasta


#remove files that are not needed
rm -f flanks.${bed_file}.${patch_reference_name}
rm -f "${first_contig}.fasta" "${second_contig}.fasta" "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa"

echo "Done."
echo "==========================="
date

