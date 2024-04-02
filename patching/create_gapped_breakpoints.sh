#!/bin/bash
#SBATCH --job-name=create_gapped_breakpoints.20240327
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=create_gapped_breakpoints.20240327.%j.log

set -e
#set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/ONT

T2T_file=$1
#T2T_file="T2T/full.maternal.contigs.T2T.scaffolds.txt"
assembly=$2
assembly_name=$(basename -- "$assembly")
assembly_name="${assembly_name%.*}"
#assembly="full.maternal.contigs.fa"
flank_size=1000000 #set the default flank length to 1Mbp
min_gap_size=10 #the minimal length of a gap worth patching

extract_flanks() {
    local coordinates="$1"
    local flank_size="$2"
    local chromosome="$3"
    local chromosome_length="$4"

    local contig
    local flank_start
    local flank_end

    contig=$(awk '{print $1}' <<< "$coordinates")
    N_start=$(awk '{print $2}' <<< "$coordinates")
    N_end=$(awk '{print $3}' <<< "$coordinates")

    #echo "N_start: $N_start"
    #echo "N_end: $N_end"

    #echo "coordinates:" $coordinates
    #echo "contig: " $contig

    flank_start=$(awk -v flank_size="$flank_size" '{
        flank_start = $2 - flank_size;
        if (flank_start < 1) {
            flank_start = 0;
        }
        print flank_start;
    }' <<< "$coordinates")

    flank_end=$(awk -v flank_size="$flank_size" -v chromosome_length="$chromosome_length" '{
        flank_end = $3 + flank_size;
        if (flank_end > chromosome_length) {
            flank_end = chromosome_length;
        }
        print flank_end;
    }' <<< "$coordinates")

    echo -e "$contig\t$flank_start\t$N_start\t$chromosome"
    echo -e "$contig\t$N_end\t$flank_end\t$chromosome"
}

#loop through all non-acrocentric contigs that contain gaps
for contig in `cat T2T/full.maternal.contigs.T2T.scaffolds.txt | tail -n +2 | egrep -v "chr13|chr14|chr15|chr21|chr22" | awk '{if ($3>0) print;}' | cut -f1`; do
    #echo $contig; #this contig contains gaps
    chromosome=`egrep $contig $T2T_file | cut -f4`
    echo $chromosome
    chromosome_length=`egrep $contig $T2T_file | cut -f2`
    
    #identify stretches of Ns
    #this is followed by merging the stretches of Ns, up to the distance of the flank (so that the flanks are always non-overlapping)
    seqtk cutN -g -n ${min_gap_size} ${assembly} | egrep ${contig} | bedtools merge -d ${flank_size} >gaps.${chromosome}.${contig}.bed

    order=1
    #for each gap, extract flanks and their coordinates
    while IFS= read -r gap || [ -n "$gap" ]; do
        if [[ -n $gap ]]; then  # Check if the line is non-empty
            bed_file="${chromosome}.gaps.${assembly_name}.${contig}.order${order}.bed"
            extract_flanks "$gap" "$flank_size" "$chromosome" "$chromosome_length" >${bed_file}
            #extract sequence as well
            bedtools getfasta -fi ${assembly} -bed ${bed_file} -name >flanks.${bed_file}.fa
            ((order++)) #increment the order if multiple gaps per chromosome
        fi
    done < gaps."$chromosome"."$contig".bed
    rm gaps.${chromosome}.${contig}.bed
done;

echo "Done."
date

