import sys
import os
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

#### code below has been copied from the following source: https://github.com/makovalab-psu/heterochromatin/blob/master/identify_short_repeats/parseTRFngsKeepHeader.py
def RotateMe(text,mode=0,steps=1):
    # function from http://www.how2code.co.uk/2014/05/how-to-rotate-the-characters-in-a-text-string/
    # Takes a text string and rotates
    # the characters by the number of steps.
    # mode=0 rotate right
    # mode=1 rotate left
    length=len(text)
     
    for step in range(steps):
    # repeat for required steps
        if mode==0:
            # rotate right
            text=text[length-1] + text[0:length-1]
        else:
            # rotate left
            text=text[1:length] + text[0]
    return text

def SmallestRotation(seq):
    smallest=seq
    for i in range(0,len(seq)):
        actual=RotateMe(seq,0,i)
        #print ("*" + actual)
        if (actual<smallest):
            #found new minimum
            smallest=actual
    return smallest

def lexicographicallySmallestRotation(seq):
    my_seq=Seq(seq)
    reverse_complement=my_seq.reverse_complement()
    reverse_complement=str(reverse_complement)

    smrt_seq=SmallestRotation(seq)
    smrt_rev_compl_seq=SmallestRotation(reverse_complement)

    #lexicographically smallest rotation is either one of the rotations of the sequence or its reverse complement
    if (smrt_seq < smrt_rev_compl_seq):
        return smrt_seq
    else:
        return smrt_rev_compl_seq

#### code above has been copied from the following source: https://github.com/makovalab-psu/heterochromatin/blob/master/identify_short_repeats/parseTRFngsKeepHeader.py

def extract_number(input_str):
    number_to_extract = re.findall(r'\d+', input_str.strip())[0] #find the first occurence
    return(int(number_to_extract))

def parse_chunks(input_str):
    chunks = input_str.strip().split('\n\n')
    parsed_chunks = []
    for chunk in chunks:
        lines = chunk.split('\n')
        parsed_chunk = {}
        new_motif= []
        for line in lines:
            key_value = line.split(' ')
            #print(len(key_value))
            if ((len(key_value) == 14) and (key_value[0].startswith("#")) and (key_value[1].startswith("position"))):

                reference = key_value[3][1]
                index = int(key_value[2]) #index of the actual motif (position within the motif)

                if (index==0): #we are starting to parse the chunk, so we can save what the motif we're parsing is
                    motif_string = previous_line.split(' ', 1)[0]
                    print(f"motif string: {motif_string}")

                mRatio = extract_number(key_value[4]) 
                threshold_mRatio=90
                #matches = extract_number(key_value[5])
                mismatches = extract_number(key_value[6])
                letterA = extract_number(key_value[9])
                letterC = extract_number(key_value[10])
                letterG = extract_number(key_value[11])
                letterT = extract_number(key_value[12])

                if ((mRatio<threshold_mRatio) and (mismatches>0)):  #if the alternative bps at certain positions are more frequent than the threshold
                    #print("mismatches are greater")
                    #print(reference)
                    print(key_value)
                    print(index)
                    #print(previous_line)
                    motif = list(motif_string) #convert into an array
                    new_motif=motif
                    #print(new_motif)
                    #print(matches)
                    #print(mismatches)

                    match reference:
                        case "A":
                             print("Reference letter is A")
                             if (letterC>letterG):
                                if (letterC>letterT):
                                    print("Replacement is letter C")
                                    new_motif[index]='C'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'
                             else:
                                if (letterG>letterT):
                                    print("Replacement is letter G")
                                    new_motif[index]='G'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'

                        case "C":
                             print("Reference letter is C")
                             if (letterA>letterG):
                                if (letterA>letterT):
                                    print("Replacement is letter A")
                                    new_motif[index]='A'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'
                             else:
                                if (letterG>letterT):
                                    print("Replacement is letter G")
                                    new_motif[index]='G'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'
                        case "G":
                             print("Reference letter is G")
                             if (letterA>letterC):
                                if (letterA>letterT):
                                    print("Replacement is letter A")
                                    new_motif[index]='A'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'
                             else:
                                if (letterC>letterT):
                                    print("Replacement is letter C")
                                    new_motif[index]='C'
                                else:
                                    print("Replacement is letter T")
                                    new_motif[index]='T'

                        case "T":
                             print("Reference letter is T")
                             if (letterA>letterC):
                                if (letterA>letterG):
                                    print("Replacement is letter A")
                                    new_motif[index]='A'
                                else:
                                    print("Replacement is letter G")
                                    new_motif[index]='G'
                             else:
                                if (letterC>letterG):
                                    print("Replacement is letter C")
                                    new_motif[index]='C'
                                else:
                                    print("Replacement is letter G")
                                    new_motif[index]='G'
                    new_motif_to_print="".join([str(i) for i in new_motif])  #convert an array back to string
                    print(f"new motif    :    {new_motif_to_print}") 
                    motif_to_be_rotated = new_motif_to_print.replace("+","").replace("-","")
                    #only write the lexicographically smallest rotation
                    rotated_motif=lexicographicallySmallestRotation(motif_to_be_rotated)
                    print(f"rotated motif:    {rotated_motif}") 
                    mo.write(rotated_motif + '\n')     #remove the extra '+' symbol at the end of the line           
            previous_line=line
            #parsed_chunk[matches] = mismatches
        #parsed_chunks.append(parsed_chunk)
    return parsed_chunks

# Open the file for reading
filename = sys.argv[1]
print(filename)
mo = open(filename + ".ncrf.new.motifs.txt", "w")

with open(filename, 'r') as f:
    # Read the entire contents of the file into a string variable
    file_contents = f.read()

parse_chunks(file_contents)
mo.close()
