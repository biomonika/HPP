#!/bin/bash
#SBATCH --job-name=extract_chromosome.20231129
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=extract_chromosome.20231129.%j.log

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

threadCount=32
minIdentity=95
reference="chm13v2.0.fa.gz"
assembly=$1
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"

mashmap --threads ${threadCount} --perc_identity ${minIdentity} --noSplit -r ${reference} -q ${assembly} -o ${assembly_name}.mashmap.tmp

cat ${assembly_name}.mashmap.tmp | awk '{if ($5=="-") print $1}' | sed 's/>//g' >${assembly_name}.contigsToBeReverseComplemented.lst
cat ${assembly_name}.mashmap.tmp | awk '{if ($5=="+") print $1}' | sed 's/>//g' >${assembly_name}.contigsToBeKeptAsIs.lst

seqtk subseq ${assembly} ${assembly_name}.contigsToBeReverseComplemented.lst >${assembly_name}.contigsToBeReverseComplemented.tmp
seqtk seq -r ${assembly_name}.contigsToBeReverseComplemented.tmp >${assembly_name}.contigsToBeReverseComplemented.fa
rm ${assembly_name}.contigsToBeReverseComplemented.tmp
seqtk subseq ${assembly} ${assembly_name}.contigsToBeKeptAsIs.lst >${assembly_name}.contigsToBeKeptAsIs.fa
        
cat ${assembly_name}.contigsToBeReverseComplemented.fa ${assembly_name}.contigsToBeKeptAsIs.fa | pigz -p ${threadCount} >${assembly_name}.unifiedAssembly.fa.gz

rm ${assembly_name}.mashmap.tmp ${assembly_name}.contigsToBeKeptAsIs.lst ${assembly_name}.contigsToBeReverseComplemented.lst ${assembly_name}.contigsToBeReverseComplemented.fa ${assembly_name}.contigsToBeKeptAsIs.fa

echo "Done."
date

