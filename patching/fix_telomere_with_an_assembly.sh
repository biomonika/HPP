#!/bin/bash
#SBATCH --job-name=fix_telomere_with_an_assembly.20240403
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --output=fix_telomere_with_an_assembly.20240403.%j.log

set -e
set -x

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
mashmap=$4

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
if [ ! -e "${mashmap}" ]; then
    echo "mashmap file does not exist. Quitting the script."
    exit 1
fi

if [ $(wc -l < "${bed_file}") -ne 1 ]; then
    echo "Bed file does not have exactly one line. Exiting script."
    exit 1
fi

echo ${bed_file} ${patch_reference} ${haplotype}

chromosome=$(echo "$bed_file" | cut -d'.' -f1)

#case "$chromosome" in
#    chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY)
#        echo "Input chromosome is valid: $chromosome"
#        ;;
#    *)
#        echo "Input chromosome is not valid. Quitting..."
#        exit 1
#        ;;
#esac

contig_to_be_patched=$(echo "$bed_file" | cut -d'.' -f1)
order=$(echo "$bed_file" | cut -d'.' -f8)

#extract flanks and find out where they belong
if [ ! -f "${bed_file}.fa" ]; then
    echo "Flank file does not exist. Creating now."
    bedtools getfasta -fi ${assembly} -bed ${bed_file} -name >${bed_file}.fa
else
    echo "Flank exists. It will not be extracted again"
fi

flank_file=${bed_file}."fa"

#BREAKPOINTS TO AN ASSEMBLY AVAILABLE FOR PATCHING

if [ -e "${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file exists and won't be re-written."
else
    echo "Wfmash does not exist. Creating now."
    wfmash --threads ${threadCount} --segment-length=1000 --map-pct-id=${minIdentity} --no-split ${patch_reference} ${flank_file} >tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    wait
    cat tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt | sed s'/\t/ /g' | cut -d' ' -f1-10 | sort -k1,1n >${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
    #remove temporary wfmash file
    #rm -f tmp.${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt
fi

if [ ! -s "${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt" ]; then
    echo "Wfmash file is empty. Flanks were not mapped. Exiting script."
    exit 1
fi

wait

file_name=("${bed_file}.${assembly_name}.TO.${patch_reference_name}.txt")

#get the missing sequence containing telomere

first_coordinate=$(awk 'NR==1 {print $8}' "${file_name}")
second_coordinate=$(awk 'NR==1 {print $9}' "${file_name}")
total_length=$(awk 'NR==1 {print $7}' "${file_name}")
contig_name_patch_reference=$(awk 'NR==1 {print $6}' "${file_name}")
contig_name_assembly=`cat ${mashmap} | egrep "\b${chromosome}\b" | cut -d' ' -f1`

distance_from_beginning=${first_coordinate}
distance_from_end=$((total_length - first_coordinate))

echo "distance_from_beginning: $distance_from_beginning"
echo "distance_from_end: $distance_from_end"

if [ "$distance_from_beginning" -lt "$distance_from_end" ]; then
    echo "we will be fixing the BEGINNING of the chromosome"
    region=${contig_name_patch_reference}:0-${first_coordinate}
    echo ${region}
    samtools faidx ${patch_reference} ${region} >${patch_reference_name}.beginning.fa
    samtools faidx ${assembly} ${contig_name_assembly} >${contig_name_assembly}.unpatched.fa
    
    #COMBINE BOTH ASSEMBLIES
    #add header
    echo ">${chromosome}.${contig_name_assembly}.telofix" >tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa
    #add sequence
    cat ${patch_reference_name}.beginning.fa ${contig_name_assembly}.unpatched.fa | seqtk seq | egrep -v "^>" | tr -d '\n' >>tmp.${chromosome}.PATCHED.${assembly_name}.telofix.fa

    #remove unnecessary files
    rm -f ${patch_reference_name}.beginning.fa ${contig_name_assembly}.unpatched.fa

else
    echo "we will be fixing the END of the chromosome"
    region=${contig_name_patch_reference}:${second_coordinate}-${total_length}
    echo ${region}
    samtools faidx ${patch_reference} ${region} >${patch_reference_name}.ending.fa
    samtools faidx ${assembly} ${contig_name_assembly} >${contig_name_assembly}.unpatched.fa

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
