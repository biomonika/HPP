import sys

if len(sys.argv) < 2:
    print("Usage: python filter_for_length.py input_file")
    sys.exit(1)

input_file = sys.argv[1]

# Generating output file name
output_file = input_file.replace('.gff', '.filtered.gff')

# Check if replacement was done
if output_file == input_file:
    raise ValueError("Input file must have a '.gff' extension")

filtered_lines = []

with open(input_file, 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
            
        # Split the line by semicolons to get key-value pairs
        pairs = line.strip().split(';')
        spacer_value = None
        
        # Extract the spacer value
        for pair in pairs:
            key, value = pair.split('=')
            if key == 'spacer':
                spacer_value = int(value)
                break
        
        # Add the line to the filtered lines list if spacer value is less than 17
        if spacer_value is None or spacer_value < 16:
            filtered_lines.append(line)

# Write the filtered lines to the output file
with open(output_file, 'w') as outfile:
    outfile.writelines(filtered_lines)

print(f"Filtered content written to {output_file}")
