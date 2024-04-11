#!/bin/bash

# This scripts replaces rDNA in the fasta file with the repeated sequence of the rDNA consensus
# input_chromosome => input chromosome that needs to have rDNA arrays replaced
# bed_file => coordinates of rDNA arrays in the input_chromosome file
# consensus_rDNA => fasta file with consensus rDNA
# experimental_length => experimentally derived length of rDNA arrays in kb

# written by Monika Cechova (mcechova@ucsc.edu)

set -e
set -x

# Function to insert rDNA into the chromosome
insert_rDNA() {
    input_chromosome="$1"
    bed_file="$2"
    consensus_rDNA="$3"
    experimental_length="$4" #set experimental length in kb

    num_lines=$(wc -l < "$bed_file")
    if [ "$num_lines" -ne 1 ]; then
        echo "Error: $bed_file must have exactly one line."
        exit 1
    fi

    # Calculate consensus rDNA length
    consensus_rDNA_length=$(bioawk -c fastx '{ print length($seq) }' < "$consensus_rDNA")

    if (( consensus_rDNA_length == 0 )); then
        echo "Error: $consensus_rDNA must have at least one sequence."
        exit 1
    fi

    # Calculate the number of copies of rDNA to insert
    number_of_copies=$(( experimental_length * 1000 / consensus_rDNA_length ))

    echo "We need to repeat consensus $number_of_copies times in order to reach the experimental length of $experimental_length kb"

    # Read the sequence from the FASTA file
    sequence=$(awk '/^>/ {next} {printf $0}' "$consensus_rDNA")

    # Repeat the sequence based on the number of repeats
    for ((i=0; i<number_of_copies; i++)); do
        repeated_sequence="${repeated_sequence}${sequence}"
    done

    # Output repeated sequence in FASTA format
    echo ">Repeated_sequence"
    echo "$repeated_sequence"

    header=$(awk 'NR==1 {sub(/^>/, ""); print; exit}' "$input_chromosome")
    input_chromosome_length=$(bioawk -c fastx '{ print length($seq) }' < "$input_chromosome")

    rDNA_start=$(awk 'NR==1 {print $2}' "$bed_file")
    rDNA_end=$(awk 'NR==1 {print $3}' "$bed_file")

    echo "header: $header"
    echo "rDNA_start: $rDNA_start"
    echo "rDNA_end: $rDNA_end"

    #going from bed to gff, increment end coordinates
    rDNA_end=$((rDNA_end+1))

    samtools faidx ${input_chromosome} ${header}:1-${rDNA_start} > first_part.${input_chromosome}
    samtools faidx ${input_chromosome} ${header}:${rDNA_end}-${input_chromosome_length} > second_part.${input_chromosome}

    echo ">"${input_chromosome}.$(echo "rDNA_padded_$header" | tr -d '>') >${input_chromosome}.padded.tmp
    cat "first_part.${input_chromosome}" | grep -v ">" | tr -d '\n' >>${input_chromosome}.padded.tmp
    echo $repeated_sequence >>${input_chromosome}.padded.tmp
    cat "second_part.${input_chromosome}" | grep -v ">" | tr -d '\n' >>${input_chromosome}.padded.tmp

    cat ${input_chromosome}.padded.tmp | seqtk seq -l 60 >${input_chromosome}.padded.fasta
    echo "rDNA padded sequence written to ${input_chromosome}.padded.fasta"

    rm -f ${input_chromosome}.padded.tmp
    rm -f first_part.${input_chromosome}
    rm -f second_part.${input_chromosome}
}

# Example usage
input_chromosome=$1
bed_file=$2
consensus_rDNA=$3
experimental_length=$4

insert_rDNA "$input_chromosome" "$bed_file" "$consensus_rDNA" "$experimental_length"
