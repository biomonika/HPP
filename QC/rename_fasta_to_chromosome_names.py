#this script takes assembly and metadata file in the form of T2Tscaffolds.txt, and renames the scaffold names to meaningful chromosomes

import os
import argparse

def load_mapping(metadata_path):
    mapping = {}
    with open(metadata_path) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                original = parts[0]
                new = parts[3]
                mapping[original] = new
    return mapping

def replace_fasta_headers(fasta_path, output_path, mapping):
    with open(fasta_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith(">"):
                parts = line[1:].strip().split(maxsplit=1)
                original_name = parts[0]
                rest = parts[1] if len(parts) > 1 else ""
                new_prefix = mapping.get(original_name, original_name)
                combined_name = f"{new_prefix}.{original_name}"
                fout.write(f">{combined_name} {rest}\n" if rest else f">{combined_name}\n")
            else:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers using scaffold-to-chromosome mapping from a metadata file."
    )
    parser.add_argument("fasta_file", help="Path to input FASTA file")
    parser.add_argument("metadata_file", help="Path to metadata file with mapping") #provide T2T.scaffolds.txt file
    args = parser.parse_args()

    base, ext = os.path.splitext(args.fasta_file)
    output_file = f"{base}.renamed{ext or '.fa'}"

    mapping = load_mapping(args.metadata_file)
    replace_fasta_headers(args.fasta_file, output_file, mapping)

    print(f"✅ Done. Output written to: {output_file}")

if __name__ == "__main__":
    main()
