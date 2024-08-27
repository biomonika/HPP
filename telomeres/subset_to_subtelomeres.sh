#!/bin/bash
#SBATCH --job-name=subset_to_subtelomeres.20240823
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=subset_to_subtelomeres.20240823.%j.log
#SBATCH --time=12:00:00

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/ONT

#the input bam file that needs to be subsetted
bam_file=$1

# Check if the input file exists
if [ ! -f "$bam_file" ]; then
    echo "Error: Input file '$bam_file' not found."
    exit 1
fi

base_name=$(basename "$bam_file" .bam)
subset_bam="${base_name}_subset.bam"

# Declare an associative array
declare -A sequence_dict

# Extract sequence names and lengths from BAM header and store in the associative array
while IFS=$' ' read -r seq_name length; do
    sequence_dict["$seq_name"]=$length
done < <(samtools view -H "$bam_file" | grep "^@SQ" | awk '{print $2, $3}' | sed 's/SN://;s/LN://')


# Create the subset BAM file by including the header first
samtools view -H "$bam_file" > "$subset_bam"


for seq_name in "${!sequence_dict[@]}"; do
    #print the sequence/chromosome name and the associated sequence length
    echo "$seq_name: ${sequence_dict[$seq_name]}"
    
    #retrieve the sequence length from the bam header
    length=${sequence_dict[$seq_name]}

    start=1
    flank=20000 #how much of the reference we want to keep -> reads overlapping this interval will be retrieved

    last_start=$((length - flank + 1))
    last_end=$length

    #we will get the coordinates at the beginning and the end of each sequence/chromosome
    echo $start "-" $flank
    echo $last_start "-" $last_end

    # Extract reads from the first and last 20 kb of each chromosome
    samtools view "$bam_file" "$seq_name:$start-$flank" >> "$subset_bam"
    samtools view "$bam_file" "$seq_name:$last_start-$last_end" >> "$subset_bam"
done

sorted_bam="${subset_bam}_sorted.bam"
samtools sort -o "$sorted_bam" "$subset_bam"
mv ${sorted_bam} ${subset_bam}
samtools index ${subset_bam}

echo ""
samtools flagstat ${subset_bam}

echo "Done."

