#!/bin/bash

set -e

source /opt/miniconda/etc/profile.d/conda.sh

assembly=$1
assembly_name="$(basename ${assembly})"
assembly_name="${assembly%.*}"

myThreads=10

conda activate /public/home/mcechova/conda/ncrf

echo "Using set of "
awk '{print NF}' StSat.txt | sort -nu | tail -n 1
echo " motifs."

#identify the motifs using the set from the heterochromatin paper as a template
motifs="$(cat "StSat.txt")"
cat ${assembly} | NCRF --stats=events --positionalevents --minlength=500 ${motifs} >NCRF.${assembly_name}.raw1.ncrf

conda activate /public/home/mcechova/conda/python.3.10
#find the candidate new motifs that should be added to the original set
python parse_ncrf_to_find_new_motifs.py NCRF.${assembly_name}.raw1.ncrf #adds the extension ".ncrf.new.motifs.txt"

#convert into one motif per line, so that only unique motifs can be kept
cat StSat.txt | tr ' ' '\n' >StSat.tmp
cat NCRF.${assembly_name}.raw1.ncrf.ncrf.new.motifs.txt StSat.tmp | sort | uniq | tr '\n' ' ' >extendedStSat.txt #concatenated original motifs and the extended motifs
rm StSat.tmp

conda activate /public/home/mcechova/conda/ncrf
echo "Using set of "
awk '{print NF}' extendedStSat.txt | sort -nu | tail -n 1
echo " motifs."
#identify the motifs using the set from the heterochromatin paper as a template
motifs="$(cat "extendedStSat.txt")"
cat ${assembly} | NCRF --stats=events --positionalevents --minlength=500 ${motifs} >NCRF.${assembly_name}.raw2.ncrf


#find the overlaps in the motifs
ncrf_cat NCRF.${assembly_name}.raw2.ncrf | ncrf_resolve_overlaps >${assembly_name}.overlaps.txt

conda activate /public/home/mcechova/conda/python.3.10
#only report the overlapping alignment with the highest mRatio
python report_highest_scoring_motif_per_overlapping_block.py ${assembly_name}.overlaps.txt
awk '{print $3 "\t" $4 "\t" $5 "\t" $2 "\t" $9 "\t" $6}' ${assembly_name}.overlaps.txt >NCRF.${assembly_name}.bed
echo "Done."
