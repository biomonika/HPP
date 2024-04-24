#!/bin/bash
#SBATCH --job-name=create_matpat_PAN027_assembly
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=7gb
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=create_matpat_PAN027_assembly.20240322%j.log
#SBATCH --time=3:00:00

#Written by Monika Cechova mcechova@ucsc.edu

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh; 
conda activate /private/home/mcechova/conda/alignment

mkdir -p PAN027/maternal
mkdir -p PAN027/paternal

### RUN USING FULL ###
hap1_assignment="parent_of_origin.parental.kmers.PAN027.haplotype1.full.verkko2.0.scaff.unifiedAssembly.txt"
hap2_assignment="parent_of_origin.parental.kmers.PAN027.haplotype2.full.verkko2.0.scaff.unifiedAssembly.txt"

#derive assembly names from the assignment file
hap1_assembly="${hap1_assignment/parent_of_origin.parental.kmers./}"
hap1_assembly="${hap1_assembly/.txt/.fa}"
hap2_assembly="${hap2_assignment/parent_of_origin.parental.kmers./}"
hap2_assembly="${hap2_assembly/.txt/.fa}"

cat ${hap1_assignment} ${hap2_assignment} | grep "maternal" | cut -f1 | sort >full.maternal.scaffolds.txt
cat ${hap1_assignment} ${hap2_assignment} | grep "paternal" | cut -f1 | sort >full.paternal.scaffolds.txt

cat ${hap1_assembly} ${hap2_assembly} | seqtk subseq - full.maternal.scaffolds.txt >PAN027/maternal/full.maternal.scaffolds.fa
cat ${hap1_assembly} ${hap2_assembly} | seqtk subseq - full.paternal.scaffolds.txt >PAN027/paternal/full.paternal.scaffolds.fa

### RUN USING DUPLEX ###
hap1_assignment="parent_of_origin.parental.kmers.PAN027.haplotype1.duplex.verkko2.0.scaff.unifiedAssembly.txt"
hap2_assignment="parent_of_origin.parental.kmers.PAN027.haplotype2.duplex.verkko2.0.scaff.unifiedAssembly.txt"

#derive assembly names from the assignment file
hap1_assembly="${hap1_assignment/parent_of_origin.parental.kmers./}"
hap1_assembly="${hap1_assembly/.txt/.fa}"
hap2_assembly="${hap2_assignment/parent_of_origin.parental.kmers./}"
hap2_assembly="${hap2_assembly/.txt/.fa}"

cat ${hap1_assignment} ${hap2_assignment} | grep "maternal" | cut -f1 | sort >duplex.maternal.scaffolds.txt
cat ${hap1_assignment} ${hap2_assignment} | grep "paternal" | cut -f1 | sort >duplex.paternal.scaffolds.txt

cat ${hap1_assembly} ${hap2_assembly} | seqtk subseq - duplex.maternal.scaffolds.txt >PAN027/maternal/duplex.maternal.scaffolds.fa
cat ${hap1_assembly} ${hap2_assembly} | seqtk subseq - duplex.paternal.scaffolds.txt >PAN027/paternal/duplex.paternal.scaffolds.fa

echo "Done."
date

