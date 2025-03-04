from Bio import SeqIO
import sys
import os

def split_fasta_by_suffix(input_fasta):
    """
    Reads a FASTA file and creates new files based on the suffix in each header.
    """
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]  # Get input filename without extension
    
    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.id
            suffix = header.split(".")[-1]  # Extract suffix after the last dot
            output_filename = f"{base_name}_{suffix}.fasta"
            
            with open(output_filename, "w") as out_fasta:
                SeqIO.write(record, out_fasta, "fasta")
            
            print(f"Written {header} to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_diploid_into_haplotypes.py input.fasta")
    else:
        split_fasta_by_suffix(sys.argv[1])
