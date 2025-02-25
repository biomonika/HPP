#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="generate_vcf_header.20250225"
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=generate_vcf_header.20250225.%j.log

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/microsatellite

ref_fasta="$1"
vcf_file="$2"

if [[ -z "$vcf_file" || -z "$ref_fasta" ]]; then
    echo "Usage: $0 <reference_fasta> <vcf_file>"
    exit 1
fi

# Extract chromosome names and lengths from the FASTA index
samtools faidx "$ref_fasta"

# Create VCF header
header="##fileformat=VCFv4.2"
while read -r line; do
    chrom=$(echo "$line" | cut -f1)
    length=$(echo "$line" | cut -f2)
    header+=$'\n'"##contig=<ID=$chrom,length=$length>"
done < "${ref_fasta}.fai"

header+=$'\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample'

# Derive the name for the compressed file
compressed_vcf="${vcf_file%.vcf}.vcf.gz"

# Combine the new header with the VCF body, sort, and compress
(echo "$header"; cat "$vcf_file" | sort -k1,1 -k2,2n) | bgzip > "$compressed_vcf"

# Index the compressed VCF
tabix -p vcf "$compressed_vcf"

# Cleanup
rm "${ref_fasta}.fai"

echo "Compressed VCF with new header saved to $compressed_vcf"
