#!/bin/bash
#SBATCH --job-name=fix_telomere_with_an_assembly.20240528
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --output=fix_telomere_with_an_assembly.20240403.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=24
minIdentity=95
bed_file=$1
assembly=$2
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"
patch_reference=$3
patch_reference_name=$(basename -- "$patch_reference")
patch_reference_name="${patch_reference_name%.*}"

#Check if input files exist

if [ ! -e "${bed_file}" ]; then
    echo "Bed file does not exist. Quitting the script."
    exit 1
fi
if [ ! -e "${assembly}" ]; then
    echo "Assembly file to be patched does not exist. Quitting the script."
    exit 1
fi
if [ ! -e "${patch_reference}" ]; then
    echo "Reference that should be used for patching does not exist. Quitting the script."
    exit 1
fi

if [ $(wc -l < "${bed_file}") -ne 1 ]; then
    echo "Bed file does not have exactly one line. Exiting script."
    exit 1
fi

echo ${bed_file} ${patch_reference} ${haplotype}

chromosome=$(echo "$bed_file" | cut -d'.' -f1)
echo "chromosome: $chromosome"

case "$chromosome" in
    chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY)
        echo "Input chromosome is valid: $chromosome"
        ;;
    *)
        echo "Input chromosome is not valid. Quitting..."
        exit 1
        ;;
esac

#chromosome name should be prefix of the bed file
#order should be included in the file name
#contig name will be extracted from the bed file
contig_to_be_patched=$(echo "$bed_file" | cut -d'.' -f1)
order=$(echo "$bed_file" | grep -oP 'order\d+')
contig_name_assembly=$(awk 'NR==1 {print $1}' "$bed_file")

#extract flanks and find out where they belong
if [ ! -f "${bed_file}.fa" ]; then
    echo "Flank file does not exist. Creating now."
    bedtools getfasta -fi ${assembly} -bed ${bed_file} -name >${bed_file}.fa
    #rename the header to only contain the contig name
    sed -i "s/^>.*/>$contig_name_assembly/" "${bed_file}.fa"
else
    echo "Flank exists. It will not be extracted again"
fi

flank_file=${bed_file}."fa"
alignment_padding_left=0 #if left flank does not align fully, this is how much of a padding there is
alignment_padding_right=0 #if right flank does not align fully, this is how much of a padding there is

#BREAKPOINTS TO AN ASSEMBLY AVAILABLE FOR PATCHING

if [ -e "${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file exists and won't be re-written."
else
    echo "Wfmash does not exist. Creating now."
    wfmash --threads ${threadCount} --segment-length=1000 --map-pct-id=${minIdentity} --no-split ${patch_reference} ${flank_file} >tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    wait
    cat tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt | sed s'/\t/ /g' | cut -d' ' -f1-10 >${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    #remove temporary wfmash file
    rm -f tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
fi

wfmash_file="${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt"

flank_size=$(awk 'NR==1 {print $2}' "$wfmash_file")
flank_start=$(awk 'NR==1 {print $3}' "$wfmash_file")
flank_end=$(awk 'NR==1 {print $4}' "$wfmash_file")

#if the flank does not align fully, the sequence should be truncated by alignment_padding_left
#this is only useful for telomeres at starts
#alignment_padding_left tells you how much bigger the patch should be
alignment_padding_left=$flank_start
#if the flank does not align fully, the sequence should be truncated by alignment_padding_right
#this is only useful for telomeres at ends
#alignment_padding_right tells you how much bigger the patch should be
alignment_padding_right=$((flank_size - flank_end)) 

echo "alignment_padding_left: $alignment_padding_left"
echo "alignment_padding_right: $alignment_padding_right"

if [ ! -s "${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file is empty. Flanks were not mapped. Exiting script."
    exit 1
fi

wait

file_name=("${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt")

#get the missing sequence containing telomere

first_coordinate=$(awk 'NR==1 {print $8}' "${file_name}") #where flank maps in reference
second_coordinate=$(awk 'NR==1 {print $9}' "${file_name}") #where flank maps in reference
#going from bed to gff, increment second_coordinate
second_coordinate=$((second_coordinate+1))

total_length=$(awk 'NR==1 {print $7}' "${file_name}")
contig_name_patch_reference=$(awk 'NR==1 {print $6}' "${file_name}")

distance_from_beginning=${first_coordinate} #distance measured using flank mapped to the reference
distance_from_end=$((total_length - second_coordinate)) #distance measured using flank mapped to the reference

echo "distance_from_beginning: $distance_from_beginning"
echo "distance_from_end: $distance_from_end"

if [[ $bed_file == *"start"* ]]; then
    echo "we will be fixing the BEGINNING of the chromosome"
    
    #flank_padding_left describes the position of the flank in relation to the assembly
    #most commonly this would be 0, meaning flank is at the edge of the assembly
    flank_coordinate_left=$(awk 'NR==1 {print $2}' "${bed_file}")
    echo "flank_coordinate_left: $flank_coordinate_left"

    echo "Adjust for padding."
    #we need to skip the distance to the flank (derived from the bed file), and also the part of the flank that is not aligned
    adjusted_flank_coordinate_left=$((flank_coordinate_left + alignment_padding_left))
    echo "adjusted_flank_coordinate_left: $adjusted_flank_coordinate_left"

    #region is extracted from the reference that contains telomeric sequence
    region=${contig_name_patch_reference}:0-${first_coordinate}
    echo "region: $region"
    echo ${patch_reference} ${region}
    echo ${assembly} ${contig_name_assembly}
    samtools faidx ${patch_reference} ${region} >${patch_reference_name}.beginning.fa

    #unpatched is the sequence from the assembly that will remain unchanged
    #if flank is not at the edge of the assembly, then we need to truncate the sequence from the assembly
    #it will be truncated by the distance between the flank and the actual assembly edge
    samtools faidx ${assembly} ${contig_name_assembly}:${adjusted_flank_coordinate_left}- >${contig_name_assembly}.unpatched.fa
    
    #COMBINE BOTH ASSEMBLIES
    #add header
    echo ">${chromosome}.${contig_name_assembly}.telofix" >tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa
    #add sequence
    cat ${patch_reference_name}.beginning.fa ${contig_name_assembly}.unpatched.fa | seqtk seq | egrep -v "^>" | tr -d '\n' >>tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa

    #remove unnecessary files
    rm -f ${patch_reference_name}.beginning.fa ${contig_name_assembly}.unpatched.fa

else
    echo "we will be fixing the END of the chromosome"

    #calculate the length of the assembly contig
    samtools faidx ${assembly} ${contig_name_assembly} >tmp.${contig_name_assembly}.unpatched.fa
    assembly_contig_length=$(bioawk -c fastx '{ print length($seq) }' < "tmp.${contig_name_assembly}.unpatched.fa")
    
    #flank_padding_right describes the position of the flank in relation to the assembly
    #most commonly this would be the length of the assembly sequence, meaning flank is at the edge of the assembly
    flank_coordinate_right=$(awk 'NR==1 {print $3}' "${bed_file}")
    adjusted_flank_coordinate_right=$((flank_coordinate_right - alignment_padding_right)) #adjust if the flank doesn't align fully

    flank_padding_right=$((assembly_contig_length - flank_coordinate_right))
    echo "flank_padding_right: $flank_padding_right"

    region=${contig_name_patch_reference}:${second_coordinate}-
    echo ${patch_reference} ${region}
    echo ${assembly} ${contig_name_assembly}
    samtools faidx ${patch_reference} ${region} >${patch_reference_name}.ending.fa

    #unpatched is the sequence from the assembly that will remain unchanged
    #if flank is not at the edge of the assembly, then we need to truncate the sequence from the assembly
    #it will be truncated by the distance between the flank and the actual assembly edge
    samtools faidx ${assembly} ${contig_name_assembly}:0-${adjusted_flank_coordinate_right} >${contig_name_assembly}.unpatched.fa
    rm tmp.${contig_name_assembly}.unpatched.fa

    #COMBINE BOTH ASSEMBLIES
    #add header
    echo ">${chromosome}.${contig_name_assembly}.telofix" >tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa
    #add sequence
    cat ${contig_name_assembly}.unpatched.fa ${patch_reference_name}.ending.fa | seqtk seq | egrep -v "^>" | tr -d '\n' >>tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa

    #remove unnecessary files
    rm -f ${patch_reference_name}.beginning.fa ${contig_name_assembly}.unpatched.fa

fi

#reformat to 60 characters per line in fasta file
seqtk seq -l 60 tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa >${chromosome}.PATCHED.${assembly_name}.telofix.fa

#remove unnecessary files
rm tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa

echo "Done."
echo "==========================="
date
