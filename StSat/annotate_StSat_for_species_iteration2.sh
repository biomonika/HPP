#!/bin/bash

set -e
set -x

source /opt/miniconda/etc/profile.d/conda.sh

assembly=$1
assembly_name="$(basename ${assembly})"
assembly_name="${assembly%.*}"

conda activate /public/home/mcechova/conda/ncrf
echo "Using set of "
awk '{print NF}' extendedStSat.txt | sort -nu | tail -n 1
echo " motifs."
#identify the motifs using the set from the heterochromatin paper as a template
motifs="$(cat "extendedStSat.txt")"
cat ${assembly} | NCRF --stats=events --positionalevents --minlength=500 ${motifs} >NCRF.${assembly_name}.raw2.ncrf

#find the overlaps in the motifs from the second iteration
ncrf_cat NCRF.${assembly_name}.raw2.ncrf | ncrf_summary --minmratio=0.85 >${assembly_name}.summary2.txt
ncrf_resolve_overlaps ${assembly_name}.summary2.txt >${assembly_name}.overlaps2.txt

conda activate /public/home/mcechova/conda/python.3.10
#only report the overlapping alignment with the highest mRatio
python report_highest_scoring_motif_per_overlapping_block.py ${assembly_name}.overlaps2.txt

awk '{print $3 "\t" $4 "\t" $5 "\t" $2 "\t" $9 "\t" $6}' ${assembly_name}.overlaps2.txt.unique.entries.txt | sort -k 1,1 -k2,2n >NCRF.${assembly_name}.iteration2.bed
cat NCRF.${assembly_name}.iteration2.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo "Done."
