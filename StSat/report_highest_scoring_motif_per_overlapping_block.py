import sys
import re

def extract_number(string_with_percentage):
    # extract percentage using regular expression
    percentage = re.findall(r'\d+\.\d+%', string_with_percentage)[0]

    # remove '%' character and convert to float
    decimal_percentage = float(percentage[:-1])/100
    return(decimal_percentage)

def parse_chunks(input_str):
    chunks = input_str.strip().split('\n\n')
    parsed_chunks = []
    for chunk in chunks:
        lines = chunk.split('\n')
        parsed_chunk = {}
        array_of_mratios = []
        for line in lines:
            key_value = line.split('\t')
            #print(len(key_value))
            if (len(key_value) == 13) and (not key_value[0].startswith("#")) :
                #print(key_value)
                motif = key_value[1]
                mRatio = extract_number(key_value[8]) 
                array_of_mratios.append({'line': line, 'mRatio': mRatio})
        max_mratio = max(array_of_mratios, key=lambda x: x['mRatio'])['mRatio'] #find the biggest mRatio in the chunk
        #print("MAX IS:")
        #print(str(max_mratio))
        for item in array_of_mratios: #loop through the dictionary to find the line associated with the biggest mRatio
            if item['mRatio'] == max_mratio:
                line_with_biggest_mratio = item['line']
                print(line_with_biggest_mratio)
                ue.write(line_with_biggest_mratio + '\n')
                break
    return parsed_chunks

# Open the file for reading
filename = sys.argv[1]
print(filename)
ue = open(filename + ".unique.entries.txt", "w")

with open(filename, 'r') as f:
    # Read the entire contents of the file into a string variable
    file_contents = f.read()

parse_chunks(file_contents)
ue.close()
