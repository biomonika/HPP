#!/bin/bash
#SBATCH --job-name=repeatmask.20230616
#SBATCH --partition=main
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --output=repeatmask.20230616.%j.log
#SBATCH --time=168:00:00

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/alignment

#LINE was identified by greping the original repeatmasker file from Gabby.
#For StSat, I used my coordinated from iteration 1 and NCRF.
#This is because for the chimpanzee, StSat was annotated as PTPCHT, so no StSat could be greped. It only found StSat in gorilla and bonobo. For consistency I decided to go with NCRF.
#I also want to mask GAP and SAT 

#notice that files in the format repeats.species.merged.bed are those with merged repeats revealing long stretches of predominantly satellites

#prepare the repeats track

#MERGE ALL REPEATS TO CREATE REPEAT TRACKS

#first convert .out format to .bed format
python RMOut-to-bed.py CHM13v2.0_ChrXY_FinalRepeatAnnotations_May2023.out chrY.human.bed
python RMOut-to-bed.py Chimpanzee_FinalRepeatAnnotations.May2023.out chrY.chimpanzee.bed
python RMOut-to-bed.py Bonobo_FinalRepeatAnnotations.May2023.out chrY.bonobo.bed
python RMOut-to-bed.py Gorilla_FinalRepeatAnnotations.May2023.out chrY.gorilla.bed
python RMOut-to-bed.py Sumatran_FinalRepeatAnnotations.May2023.out chrY.sorang.bed
python RMOut-to-bed.py Bornean_FinalRepeatAnnotations.May2023.out chrY.borang.bed
python RMOut-to-bed.py Siamang_FinalRepeatAnnotations.May2023.out chrY.gibbon.bed

#only large consecutive regions larger than 250000 will be kept and labeled as repeats
names=("human" "chimpanzee" "bonobo" "gorilla" "sorang" "borang" "gibbon")
for i in {1..7}; do echo $i; chromosome=${chromosomes[i]}; cat chrY.${names[i]}.bed | bedtools sort | bedtools merge -d 1000 | awk '{ if (($3-$2)>250000) print;}' | sort -V -k1,1 -k2,2 >repeats.${names[i]}.merged.bed; done;

#prepare the track for circos as well
cat chrY.human.bed | sed -e 's/chrY/hs1/' >merged.repeats.species.bed.tmp
cat chrY.chimpanzee.bed | sed -e 's/chrY/hs2/' >>merged.repeats.species.bed.tmp
cat chrY.bonobo.bed | sed -e 's/chrY/hs3/' >>merged.repeats.species.bed.tmp
cat chrY.gorilla.bed | sed -e 's/chrY/hs4/' >>merged.repeats.species.bed.tmp
cat chrY.sorang.bed | sed -e 's/chrY/hs5/' >>merged.repeats.species.bed.tmp
cat chrY.borang.bed | sed -e 's/chrY/hs6/' >>merged.repeats.species.bed.tmp
cat chrY.gibbon.bed | sed -e 's/chrY/hs7/' >>merged.repeats.species.bed.tmp

cat merged.repeats.species.bed.tmp | awk '{print $0 "\tfill_color=purple"}' | sed 's/ /\t/g' | sort -V -k1,1 -k2,2 >merged.repeats.species.bed
rm merged.repeats.species.bed.tmp


#GATHER THE HETEROCHROMATIC, REPETITVE PORTIONS OF THE Y, SO THAT THEY CAN BE MASKED AND ONLY EUCHROMATIN LEFT
cat chrY.human.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.human.bed
cat chrY.chimpanzee.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.chimpanzee.bed
cat chrY.bonobo.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.bonobo.bed
cat chrY.gorilla.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.gorilla.bed
cat chrY.sorang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.sorang.bed
cat chrY.borang.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.borang.bed
cat chrY.gibbon.bed | egrep --ignore-case --color "SAT|GAP|LINE" >sat_gap_interspersed.chrY.gibbon.bed


#WIDE MASKING, AFTER WHICH ONLY EUCHROMATIC PROPORTION SHOULD REMAIN
cat repeats.human.merged.bed sat_gap_interspersed.chrY.human.bed HSAT.human.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >human.WideHardmasked.bed
cat repeats.chimpanzee.merged.bed sat_gap_interspersed.chrY.chimpanzee.bed StSat.chimpanzee.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >chimpanzee.WideHardmasked.bed
cat repeats.bonobo.merged.bed sat_gap_interspersed.chrY.bonobo.bed StSat.bonobo.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >bonobo.WideHardmasked.bed
cat repeats.gorilla.merged.bed sat_gap_interspersed.chrY.gorilla.bed StSat.gorilla.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gorilla.WideHardmasked.bed
cat repeats.sorang.merged.bed sat_gap_interspersed.chrY.sorang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >sorang.WideHardmasked.bed
cat repeats.borang.merged.bed sat_gap_interspersed.chrY.borang.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >borang.WideHardmasked.bed
cat repeats.gibbon.merged.bed sat_gap_interspersed.chrY.gibbon.bed | cut -f1-3 | sort -k 1,1 -k2,2n | bedtools merge >gibbon.WideHardmasked.bed

#hardmask the problematic repeat regions
bedtools maskfasta -fi chrY.human.fasta -bed human.WideHardmasked.bed -fo chrY.human.WideHardmasked.fasta
bedtools maskfasta -fi chrY.chimpanzee.fasta -bed chimpanzee.WideHardmasked.bed -fo chrY.chimpanzee.WideHardmasked.fasta
bedtools maskfasta -fi chrY.bonobo.fasta -bed bonobo.WideHardmasked.bed -fo chrY.bonobo.WideHardmasked.fasta
bedtools maskfasta -fi chrY.gorilla.fasta -bed gorilla.WideHardmasked.bed -fo chrY.gorilla.WideHardmasked.fasta
bedtools maskfasta -fi chrY.sorang.fasta -bed sorang.WideHardmasked.bed -fo chrY.sorang.WideHardmasked.fasta
bedtools maskfasta -fi chrY.borang.fasta -bed borang.WideHardmasked.bed -fo chrY.borang.WideHardmasked.fasta
bedtools maskfasta -fi chrY.gibbon.fasta -bed gibbon.WideHardmasked.bed -fo chrY.gibbon.WideHardmasked.fasta


echo "Done."
