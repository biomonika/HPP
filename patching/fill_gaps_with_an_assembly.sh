#!/bin/bash
#SBATCH --job-name=fill_gaps_with_an_assembly.sh.20240329
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --output=fill_gaps_with_an_assembly.sh.20240329.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=24
minIdentity=95
reference="chm13v2.0.fa.gz"
flank_file=$1
flank_name=$(basename -- "$flank_file")
flank_name="${flank_name%.*}"
assembly=$2
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"
patch_reference=$3
patch_reference_name=$(basename -- "$patch_reference")
patch_reference_name="${patch_reference_name%.*}"
bed_file=$4
flank_size=1000000

echo ${flank_name} ${patch_reference} ${haplotype}

chromosome=$(echo "$flank_file" | cut -d'.' -f1)
order=$(echo "$flank_file" | cut -d'.' -f8)

#find out where flanks belong


#MAP BREAKPOINTS TO AN ASSEMBLY AVAILABLE FOR PATCHING

if [ -e "${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file exists and won't be re-written."
else
    echo "Wfmash does not exist. Creating now."
    wfmash --threads ${threadCount} --segment-length=1000 --map-pct-id=${minIdentity} --no-split ${patch_reference} ${flank_file} >tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
    cat tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt | sed s'/\t/ /g' | cut -d' ' -f1-10 | sort -k1,1n >${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
    #remove temporary wfmash file
    rm #tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt

    #wfmash is replacing previous version with approximate mapping and mashmap
    #mashmap --filter_mode one-to-one --threads ${threadCount} --perc_identity ${minIdentity} --noSplit --segLength 1000 -r ${patch_reference} -q flanks.${flank_name}.${assembly_name} -o ${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt > /dev/null 2>&1
fi

wait

file_names=("${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt")

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
    rm -f flanks.${flank_name}.${assembly_name}
    exit -1
fi

echo ""

#Get coordinates of the GAP/BREAKPOINT REGION
gap_start=""
gap_end=""
contig_name=""
extract_variables() {
    #This code parses wfmash output
    file_name=$1
    #FIRST FLANK
    # Extract 9th column from the first row
    gap_start=$(awk 'NR==1 {print $9}' "$file_name")
    flank_end=$(awk 'NR==1 {print $4}' "$file_name")
    adjustment=$(($flank_size - $flank_end))

    if (( adjustment > 0 )); then
        gap_start=$(($gap_start - $adjustment))
        echo "Adjusted gap_start: $gap_start"
    else
        echo "No adjustment needed. The first flank maps fully at the end."
    fi

    # Extract 10th column from the second row
    gap_end=$(awk 'NR==2 {print $8}' "$file_name")
    flank_start=$(awk 'NR==2 {print $3}' "$file_name")
    adjustment=$flank_start

    if (( adjustment > 0 )); then
        gap_end=$(($gap_end - $adjustment))
        echo "Adjusted gap_end: $gap_end"
    else
        echo "No adjustment needed. The second flank maps fully at the beginning."
    fi
    
    # Extract contig name
    # Make sure the flanks both map to the same sequence
    contig_name_row1=$(awk 'NR==1 {print $6}' "$file_name")
    contig_name_row2=$(awk 'NR==2 {print $6}' "$file_name")

    if [ "$contig_name_row1" != "$contig_name_row2" ]; then
    echo "The flanks map to two different contigs, thus their mapping cannot be used for patching. Exiting script."
    exit 1
    fi

    contig_name=$contig_name_row1

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
    echo "Negative gap sizes are currently not implemented. Exiting the script."
    exit $gap_size
    #echo "The gap_size is negative. We will be extracting reverse complement of the sequence."
    #region=${contig_name}:${gap_end}-${gap_start}
    #samtools faidx ${patch_reference} ${region} >to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
    #seqtk seq -r to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa >${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
    #rm to_be_reversed.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
else
    echo "The gap_size is positive."
    region=${contig_name}:${gap_start}-${gap_end}
    samtools faidx ${patch_reference} ${region} >${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
fi

# Extract sequence names of the contigs we will be merging
echo "Extract sequence names of the contigs we will be merging"

#this is the contig with gaps that needs gap filling
original_contig=$(awk 'NR==1 {print $1}' "${bed_file}")

samtools faidx ${assembly} ${original_contig}:0-${gap_start} > first_part.${original_contig}.fasta
samtools faidx ${assembly} ${original_contig}:${gap_end}- > second_part.${original_contig}.fasta

# Use head to extract the first line of each file
header1=${original_contig}:0-${gap_start}
header2=$(head -n 1 "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa")
header3=${original_contig}:${gap_end}-

#MERGE SEPARATELY HEADER AND THE BODY/SEQUENCE
echo ">"${chromosome}.$(echo "original_$header1" | tr -d '>')"+"$(echo "${patch_reference_name}.${gap_size}.${header2}" | tr -d '>')"+"$(echo "original_$header3" | tr -d '>') >${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "first_part.${original_contig}.fasta" | grep -v ">" >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa" | grep -v ">" >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "second_part.${original_contig}.fasta" | grep -v ">" >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp

#reformat the final patched fasta chromosome
seqtk seq -L 0 ${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp >PATCHED.${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta
rm -f ${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
echo "Merged concatenated fasta for ${chromosome} saved to" PATCHED.${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta


#remove files that are not needed
#rm -f flanks.${flank_file}.${assembly_name}
rm -f "first_part.${original_contig}.fasta" "second_part.${original_contig}.fasta" "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa"

echo "Done."
echo "==========================="
date

