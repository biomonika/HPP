import sys
from itertools import combinations

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 2:
    print("Usage: python script.py input_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = "output_" + input_file

data_dict = {}

with open(input_file, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')  # Assuming columns are tab-separated; adjust accordingly

        key = columns[0]
        values = columns[1:]

        #use read name as the key
        if key in data_dict:
            data_dict[key].append(values)
        else:
            data_dict[key] = [values]

with open(output_file, 'w') as output:
    for key, values in data_dict.items():
        if len(values) >= 2:
           #if the same read name is present at least twice, create all pairwise combinations 
           pairs = list(combinations(values, 2))
           for pair in pairs:
                formatted_pair = "\t".join(pair[0] + pair[1])
                #print(f"Pairwise combination for {key}: {pair}")
                output.write(f"{formatted_pair}\n")
