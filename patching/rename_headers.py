import os
import sys

def add_newline_if_missing(file_path):
    with open(file_path, 'rb+') as file:
        file.seek(0, os.SEEK_END)
        if file.tell() > 0:
            file.seek(-1, os.SEEK_END)
            last_char = file.read(1)
            if last_char != b'\n':
                file.write(b'\n')

def rename_fasta_header(file_path):
    # Add a newline at the end if it doesn't exist
    add_newline_if_missing(file_path)

    # Get the filename without extension
    file_name = os.path.basename(file_path)
    header = os.path.splitext(file_name)[0]

    # Read the original FASTA file
    with open(file_path, 'r') as original_file:
        fasta_content = original_file.readlines()

    # Write the updated header and original content back to the same file
    with open(file_path, 'w') as updated_file:
        updated_file.write(f">{header}\n")
        updated_file.writelines(fasta_content[1:])

    print(f"Header renamed and saved to {file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rename_headers.py path/to/your/file.fasta")
        sys.exit(1)

    file_path = sys.argv[1]
    rename_fasta_header(file_path)
