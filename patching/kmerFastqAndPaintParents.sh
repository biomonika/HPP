#!/bin/bash
#SBATCH --job-name=kmerFastqAndPaintParents
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=70gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=kmerFastqAndPaintParents.20240320%j.log
#SBATCH --time=3:00:00

#This script calculates the parent-specific k-mers and later classifies the assembly into maternal and paternal contigs using parental k-mers from PAN010 and PAN011
#Written by Monika Cechova mcechova@ucsc.edu

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh; 
conda activate /private/home/mcechova/conda/alignment

assembly=$1
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"

## create the count DB from the parents

if test -d "maternal.k21.meryl"; then
    echo "Maternal folder exists."
else
    echo "Maternal folder does not exist."
    meryl \
    threads=32 \
    count \
    compress \
    k=21 \
    PAN010.cell.HTMKJDSX2_TCGCTCTGGA-TACTCAGCTG_L004_R1.fastq.gz \
    output maternal.k21.meryl
fi

if test -d "paternal.k21.meryl"; then
    echo "Paternal folder exists."
else
    echo "Paternal folder does not exist."
    meryl \
    threads=32 \
    count \
    compress \
    k=21 \
    PAN011.cell.HTMKJDSX2_CGGCGTTCTA-GAGTTGAGGC_L004_R1.fastq.gz \
    output paternal.k21.meryl
fi

## Filter for unique kmers in each dataset

if test -d "maternal.unique.meryl"; then
    echo "Maternal unique folder exists."
else
    echo "Maternal unique folder does not exist."
    meryl print \
    difference \
    maternal.k21.meryl paternal.k21.meryl.${assembly_name} \
    output maternal.unique.meryl \
    > maternal.k21.unique.dump
fi

if test -d "paternal.unique.meryl"; then
    echo "Paternal unique folder exists."
else
    echo "Paternal unique folder does not exist."
    meryl print \
    difference \
    paternal.k21.meryl maternal.k21.meryl.${assembly_name} \
    output paternal.unique.meryl \
    > paternal.k21.unique.dump
fi


## Now that we have unique parental k-mers, let's localize them in the assemblies
## Each haplotype and assembly will produce a different bed file 
## Lookup unique k-mers and create a bed file w/ their locations
meryl-lookup \
    -sequence ${assembly} \
    -mers maternal.unique.meryl \
    -output maternal.unique.meryl.${assembly_name}.bed \
    -bed \
    > maternal.unique.meryl.dump2.${assembly_name}

meryl-lookup \
    -sequence ${assembly} \
    -mers paternal.unique.meryl \
    -output paternal.unique.meryl.${assembly_name}.bed \
    -bed \
    > paternal.unique.meryl.dump2.${assembly_name}

#at this point we have unique maternal and paternal k-mers and their locations

#let's count how many such unique maternal and paternal k-mers are there per each contig
bedtools merge -i maternal.unique.meryl.${assembly_name}.bed | awk '{ lengths[$1] += $3 - $2 } END { for (contig in lengths) print contig, lengths[contig] }' | sort >maternal.kmers.${assembly_name}.txt
bedtools merge -i paternal.unique.meryl.${assembly_name}.bed | awk '{ lengths[$1] += $3 - $2 } END { for (contig in lengths) print contig, lengths[contig] }' | sort >paternal.kmers.${assembly_name}.txt

#assign contigs into maternal or paternal contigs
Rscript classify_maternal_paternal.R maternal.kmers.${assembly_name}.txt paternal.kmers.${assembly_name}.txt

#remove unnecessary files
rm -f maternal.k21.unique.dump.${assembly_name} paternal.k21.unique.dump.${assembly_name} maternal.unique.meryl.dump2.${assembly_name} paternal.unique.meryl.dump2.${assembly_name} 

#remove bed files with the location of parent-specific k-mers
#rm -f maternal.unique.meryl.${assembly_name}.bed paternal.unique.meryl.${assembly_name}.bed

#remove k-mer files
#rm -f maternal.kmers.${assembly_name}.txt paternal.kmers.${assembly_name}.txt

echo "Done."
date

