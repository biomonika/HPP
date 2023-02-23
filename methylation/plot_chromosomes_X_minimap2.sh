#!/bin/bash

set -x 
set -e

source /opt/miniconda/etc/profile.d/conda.sh; 
conda activate /public/home/mcechova/conda/modbam


#BLUE
modbamtools plot chrX.grandfather.blue.minimap2.bam -r haplotype2-0000594:58500000-60000000 --out . --prefix grandfather.blue.minimap2 --samples PAN011 &
modbamtools plot chrX.mother.blue.minimap2.bam -r haplotype2-0000887:58500000-60000000 --out . --prefix mother.blue.minimap2 --samples PAN027 &

#compare with corrected BLUE data from secphase
modbamtools plot chrX.grandfather.blue.minimap2.corrected.bam -r haplotype2-0000594:58500000-60000000 --out . --prefix grandfather.blue.minimap2.corrected --samples PAN011 &
modbamtools plot chrX.mother.blue.minimap2.corrected.bam -r haplotype2-0000887:58000000-62000000 --out . --prefix mother.blue.minimap2.corrected --samples PAN027 &

#RED
modbamtools plot chrX.grandmother.red.minimap2.bam -r haplotype1-0000042:0-4000000 --out . --prefix grandmother.red.minimap2 --samples PAN010 &
modbamtools plot chrX.mother.red.minimap2.bam -r haplotype1-0000017:57000000-61000000 --out . --prefix mother.red.minimap2 --samples PAN027 &
modbamtools plot chrX.granddaughter.red.minimap2.bam -r haplotype1-0000015:57000000-61000000 --out . --prefix granddaughter.red.minimap2 --samples PAN028 &

wait
echo "Done."
