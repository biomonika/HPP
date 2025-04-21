#!/bin/bash

#parse .dat files from TRF and create bed files with coordinates of mono- and di- nucleotides

for file in *.dat; do
    # Extract the prefix before ".part"
    prefix=$(echo "$file" | sed -E 's/(.*)\.part.*/\1/')

    # Define output file names
    mono_out="${prefix}.mono.trf.bed"
    di_out="${prefix}.di.trf.bed"

    # Extract sequence name from line 9
    sequence=$(sed -n '9p' "$file" | awk '{print $2}')

    # Process file, keeping rows from line 16 onwards
    awk -v seq="$sequence" -v mono="$mono_out" -v di="$di_out" '
        NR >= 16 {
            type = (length($14) == 1) ? mono : di;
            print seq "\t" $1-1 "\t" $2 "\t" $3 * $4 "\t" $14 >> type;
        }
    ' "$file"

    echo "Processed: $file → $mono_out, $di_out"
done

echo "Processing complete."
