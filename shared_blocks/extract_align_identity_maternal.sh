#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="extract_align_identity.20241023"
#SBATCH --cpus-per-task=32
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --partition=long
#SBATCH --output=extract_align_identity.20241023.%j.log

set -e
#set -x

numCPU=3

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

# Input files and references
tsv_file="shared_contigs_maternal_coordinates.tsv"
ref_PAN010="assembly.v1.0.PAN010.diploid.fa"
ref_pan027="assembly.v1.0.PAN027.diploid.fa"
ref_pan028="assembly.v1.0.PAN028.diploid.fa"

# Output directory for extracted sequences and alignment
output_dir="maternal_blocks"
mkdir -p "$output_dir"

# Temporary file for sequences from each haplotype
pan010_seq="$output_dir/PAN010.fa"
pan027_seq="$output_dir/PAN027.fa"
pan028_seq="$output_dir/PAN028.fa"

# Loop through each line of the TSV file
while IFS= read -r line; do

  echo "==================="
  echo $line
  echo "==================="

  tab_count=$(grep -o $'\t' <<< "$line" | wc -l)
  # Skip the row if the number of fields is not exactly 4
  if [[ $tab_count -ne 3 ]]; then
    echo "Skipping line due to incorrect number of fields"
    continue
  fi

  IFS=$'\t' read -r contig coord1 coord2 coord3 <<< "$line"

  # Parse coordinates from each haplotype (PAN010, PAN027, PAN028)
  # Assuming the format is `chromosome:start-end`
  hap1_chr=$(echo "$coord1" | cut -d ':' -f 1)
  hap1_range=$(echo "$coord1" | cut -d ':' -f 2)

  hap2_chr=$(echo "$coord2" | cut -d ':' -f 1)
  hap2_range=$(echo "$coord2" | cut -d ':' -f 2)

  hap3_chr=$(echo "$coord3" | cut -d ':' -f 1)
  hap3_range=$(echo "$coord3" | cut -d ':' -f 2)

  # Extract sequences from each haplotype using samtools faidx
  samtools faidx "$ref_PAN010" "$hap1_chr:$hap1_range" > "$pan010_seq"
  samtools faidx "$ref_pan027" "$hap2_chr:$hap2_range" > "$pan027_seq"
  samtools faidx "$ref_pan028" "$hap3_chr:$hap3_range" > "$pan028_seq"

  # Align the sequences using minimap2
  alignment_file="$output_dir/${contig}_generation1_aligned.paf"
  stat_file="$output_dir/${contig}_generation1_stat.txt"
  minimap2 -t ${numCPU} --cs -x asm5 -c "$pan027_seq" "$pan010_seq" > "$alignment_file" 2> /dev/null
  paftools.js stat "$alignment_file" >${stat_file}
  #rm $alignment_file

  alignment_file="$output_dir/${contig}_generation2_aligned.paf"
  stat_file="$output_dir/${contig}_generation2_stat.txt"
  minimap2 -t ${numCPU} --cs -x asm5 -c "$pan027_seq" "$pan028_seq" > "$alignment_file" 2> /dev/null
  paftools.js stat "$alignment_file" >${stat_file}
  #rm $alignment_file

  alignment_file="$output_dir/${contig}_generations13_aligned.paf"
  stat_file="$output_dir/${contig}_generations13_stat.txt"
  minimap2 -t ${numCPU} --cs -x asm5 -c "$pan010_seq" "$pan028_seq" > "$alignment_file" 2> /dev/null
  paftools.js stat "$alignment_file" >${stat_file}
  #rm $alignment_file

  #delete the unnecessary files
  rm "$pan010_seq" "$pan027_seq" "$pan028_seq"

done < "$tsv_file"

echo "Processing completed. Results stored in $output_dir"
