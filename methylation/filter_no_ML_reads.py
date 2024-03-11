# This script separates reads to two bam files: those that contain and do not contain Ml tag
# 2024 Monika Cechova mcechova@ucsc.edu

import sys
import pysam

filename = sys.argv[1]
samfile = pysam.AlignmentFile(filename, "rb")
output_filtered = pysam.AlignmentFile(filename+".onlyMlTags.bam", "wb", template=samfile)
output_missing = pysam.AlignmentFile(filename+".missingTags.bam", "wb", template=samfile)

counter=0
for read in samfile:
    #print(read.query_name)
    try:
        ml_tag = read.get_tag("Ml")
        output_filtered.write(read)
    except KeyError:
        #print("The 'Ml' tag was not found.")
        output_missing.write(read)
        counter+=1


samfile.close()
print("Done. Output file created. The number of reads with missing Ml tags: ", counter)
print("Please modify the script if you wish to filter ML and not Ml tags. ")
