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
bed_file="${flank_file%.fa}" #the original bed file can be derived from flank file name
flank_size=1000000

echo ${flank_name} ${patch_reference} ${haplotype}

chromosome=$(echo "$bed_file" | cut -d'.' -f1)
order=$(echo "$bed_file" | grep -o 'order[0-9]\+' | sed 's/order//')
echo "chromosome: $chromosome"
echo "order: $order"

adjust_coordinates_and_try_again () {
    #if gap filling with current coordinates is not successfull, this function should be called
    #read information about the first flank
    chromosome_of_first_flank=$(awk 'NR==1 {print $1}' "${bed_file}")
    first_coordinate_of_first_flank=$(awk 'NR==1 {print $2}' "${bed_file}")
    second_coordinate_of_first_flank=$(awk 'NR==1 {print $3}' "${bed_file}")
    name_of_first_flank=$(awk 'NR==1 {print $4}' "${bed_file}")

    #read information about the second flank
    chromosome_of_second_flank=$(awk 'NR==2 {print $1}' "${bed_file}")
    first_coordinate_of_second_flank=$(awk 'NR==2 {print $2}' "${bed_file}")
    second_coordinate_of_second_flank=$(awk 'NR==2 {print $3}' "${bed_file}")
    name_of_second_flank=$(awk 'NR==2 {print $4}' "${bed_file}")

    #adjust the coordinates, subtract for the first flank
    first_coordinate_of_first_flank=$(($first_coordinate_of_first_flank - flank_size))
    second_coordinate_of_first_flank=$(($second_coordinate_of_first_flank - flank_size))

    #adjust the coordinates, add to the seconf flank
    first_coordinate_of_second_flank=$(($first_coordinate_of_second_flank + flank_size))
    second_coordinate_of_second_flank=$(($second_coordinate_of_second_flank + flank_size))

    first_row_of_bed_file="${chromosome_of_first_flank}\t${first_coordinate_of_first_flank}\t${second_coordinate_of_first_flank}\t${name_of_first_flank}"
    second_row_of_bed_file="${chromosome_of_second_flank}\t${first_coordinate_of_second_flank}\t${second_coordinate_of_second_flank}\t${name_of_second_flank}"

    new_order=$((order+1000))
    new_bed_file_name=`echo $bed_file | sed "s/order${order}/order${new_order}/"`
    echo "new_bed_file_name: ${new_bed_file_name}"

    if [[ ${new_order} -gt 5000 ]] ; then
        echo "Giving up on gap filling. Tried too many times: ${new_order}"
        exit 1
    else
        echo "Let's try gap filling again with the new adjusted coordinates."

        echo -e ${first_row_of_bed_file} >${new_bed_file_name}
        echo -e ${second_row_of_bed_file} >>${new_bed_file_name}

        #get new flanks based on these new coordinates
        bedtools getfasta -fi ${assembly} -bed ${new_bed_file_name} >${new_bed_file_name}.fa
        sleep 5

        #calculate the average length of the extracted flank
        average_length_of_flank=`bioawk -c fastx '{ print $name, length($seq) }' ${new_bed_file_name}.fa | awk '{ total += $2 } END { print total/NR }'`
        
        if [ $average_length_of_flank -eq $flank_size ]; then
            echo "sbatch fill_gaps_with_an_assembly.sh ${new_bed_file_name}.fa ${assembly} ${patch_reference}"
            
            #remove wfmash file as it won't be needed, as it doesn't work for the match
            rm -f ${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt

            #let's re-run the script with new parameters and try again
            sbatch fill_gaps_with_an_assembly.sh ${new_bed_file_name}.fa ${assembly} ${patch_reference}
            exit 1
        else
            echo "The extracted flanks have unexpected lengths. Giving up on gap filling."
            echo "flank_size: $flank_size"
            echo "average_length_of_flank: $average_length_of_flank"
            exit 1
        fi
fi
}


#find out where flanks belong

#MAP BREAKPOINTS TO AN ASSEMBLY AVAILABLE FOR PATCHING

if [ -e "${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file exists and won't be re-written."
else
    echo "Wfmash does not exist. Creating now."
    wfmash --threads ${threadCount} --segment-length=1000 --map-pct-id=${minIdentity} --no-split ${patch_reference} ${flank_file} >tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
    cat tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt | sed s'/\t/ /g' | cut -d' ' -f1-10 | sort -k1,1n >${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
    #remove temporary wfmash file
    rm tmp.${flank_name}.${assembly_name}.TO.${patch_reference_name}.txt
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
        echo "==========================="
        echo "One of the flanks is not mapped successfully."
        #call bash function to find out if adjusting the coordinates could help
        adjust_coordinates_and_try_again
        exit 13
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
    echo "\n"
    echo "GAP/BREAKPOINT contig name (6th column): $contig_name"
    echo "GAP/BREAKPOINT start (9th column from the first row): $gap_start"
    echo "GAP/BREAKPOINT end (8th column from the second row): $gap_end"
    echo "\n"
}

# Extract coordinates from the reference used for patching
extract_variables "$max_file"

gap_size=$((gap_end - gap_start))
echo "Gap size is $gap_size"

region=""
# Check if the gap size is negative
if [ "$gap_size" -lt 0 ]; then
    echo "==========================="
    echo "Negative gap sizes ($gap_size) are currently not implemented. Exiting the script."
    
    #call bash function to find out if adjusting the coordinates could help
    adjust_coordinates_and_try_again
    exit 13
else
    echo "The gap_size is positive."
    region=${contig_name}:${gap_start}-${gap_end}
    samtools faidx ${patch_reference} ${region} >${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa
    #check to make sure the patch doesn't contain Ns
    gap_file="gaps_in_patch.${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa"
    seqtk cutN -g -n 10 ${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa >${gap_file}
    
    if [ ! -s "${gap_file}" ]; then
        echo "No gaps found in the patch."
        rm -f ${gap_file}
    else
        echo "The number of gaps inside the patch is not empty. Exiting..."
        rm -f ${gap_file}
        exit 1
    fi
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
cat "first_part.${original_contig}.fasta" | grep -v ">" | tr -d '\n' >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa" | grep -v ">" >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
cat "second_part.${original_contig}.fasta" | grep -v ">" | tr -d '\n' >>${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp

#reformat the final patched fasta chromosome
cat ${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp | seqtk seq -l 60 >${chromosome}.PATCHED.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta
rm -f ${chromosome}.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta.tmp
echo -e "\n"
echo "Merged concatenated fasta for ${chromosome} saved to" ${chromosome}.PATCHED.${order}.patched.${assembly_name}.with.${patch_reference_name}.fasta
echo "PATCHING SUCCESSFULL."

#remove files that are not needed
#rm -f flanks.${flank_file}.${assembly_name}
rm -f "first_part.${original_contig}.fasta" "second_part.${original_contig}.fasta" "${chromosome}.${order}.patch.${patch_reference_name}.${region}.fa"

echo "Done."
echo "==========================="
date

