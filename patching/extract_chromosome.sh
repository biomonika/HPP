#!/bin/bash
#SBATCH --job-name=extract_chromosome.20240329
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=extract_chromosome.20240329.%j.log

set -e
#set -x

pwd; hostname; date

# Check if the number of arguments is not equal to 2
if [ "$#" -ne 2 ]; then
    echo "Usage: extract_chromosome.sh assembly scaffolds"
    exit 1  # Exit with a non-zero status indicating an error
fi

#activate the conda environment where mashmap is installed
source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=32
minIdentity=95 #we also want tricky small contigs to map
reference="chm13v2.0.fa.gz"
assembly=$1
scaffolds=$2

assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"

mkdir -p output.${assembly_name}
echo ${assembly} ${scaffolds}

if [ -e "${assembly_name}.mashmap.txt" ]; then
    echo "File exists: $filename"
else
    echo "File does not exist: $filename"
    mashmap --threads ${threadCount} --perc_identity ${minIdentity} --segLength 1000 --noSplit -r ${reference} -q ${assembly} -o ${assembly_name}.mashmap.txt.tmp
    cat ${assembly_name}.mashmap.txt.tmp | sort -k6,6 -V -s >${assembly_name}.mashmap.txt
    rm ${assembly_name}.mashmap.txt.tmp
fi

# Array of chromosomes
chromosomes=( {1..22} X Y )

for chrom in "${chromosomes[@]}"; do
    chromosome="chr${chrom}"
    echo "Processing chromosome $chromosome"

    if egrep -q "\\b${chromosome}\\b" ${scaffolds}; then
        echo ${chromosome} " is a complete chromosome."
        contig=`egrep "\\b${chromosome}\\b" ${scaffolds} | cut -f1`
        echo ${contig}
        samtools faidx ${assembly} ${contig} >output.${assembly_name}/${chromosome}.complete.fa
    else
        echo "Not a complete chromosome. Let's check if it's truncated or if it can be patched."

        num_entries=$(egrep "\\b${chromosome}\\b" ${assembly_name}.mashmap.txt | wc -l)
        # Check if at least one file exists
        if [ "$num_entries" -gt 1 ]; then
            echo "Chromosome should be patched."
            egrep "\\b${chromosome}\\b" ${assembly_name}.mashmap.txt >${chromosome}.${assembly_name}.mashmap.txt
            Rscript --vanilla --slave parse_mashmap.R ${chromosome}.${assembly_name}.mashmap.txt > R_output.txt 2>&1
            wait
            rm ${chromosome}.${assembly_name}.mashmap.txt
        else
            echo "No files found for $chromosome"
            echo ${chromosome} >output.${assembly_name}/${chromosome}.truncated.fa
        fi
    fi
done

#sequences defined in bed files should also be writen to fasta files
for bed in chr*.${assembly_name}.mashmap.txt*bed; do 
    echo $bed
    bedtools getfasta -fi ${assembly} -bed ${bed} -name >flanks.${bed}.fa
done

echo "Done."
date
