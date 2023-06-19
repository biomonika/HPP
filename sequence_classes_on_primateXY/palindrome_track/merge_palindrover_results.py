#This script converts the default output of palindrover with individual arms into new annotation where both arms and spacer are merged into a single entry

import sys

if len(sys.argv) > 1:
    input_name = sys.argv[1]
    print("Input name:", input_name)
    with open(input_name, 'r') as input_file, open((input_name+".palindrover.bed"), 'w') as output_file:
        lines = input_file.readlines()
        for i in range(0, len(lines), 2):
            #print(i)
            line1 = lines[i].strip().split('\t')
            line2 = lines[i+1].strip().split('\t')
            #print(line1)
            #print(line2)
            merged_line = (line1[0] + "\t" + line1[1] + "\t" + line2[2])
            #print(merged_line)
            output_file.write(merged_line + '\n')
    print("Output name:", (input_name+".palindrover.bed"))
else:
    print("No input name provided.")
