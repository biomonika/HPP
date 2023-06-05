import sys
import pysam

filename = sys.argv[1]
samfile = pysam.AlignmentFile(filename, "rb")
output_5hmC = pysam.AlignmentFile(filename+".output_5hmC.bam", "wb", template=samfile)
output_5mC = pysam.AlignmentFile(filename+".output_5mC.bam", "wb", template=samfile)

for read in samfile:
    #print(read.query_name)
    mm_tags=read.get_tag("MM")
    #print(mm_tags)

    #split 5hmC and 5mC, using semicolon as a separator
    split_mm=mm_tags.split(';')
    MM_h=split_mm[0] #contains 5hmC mm tag, positions of methylations (how many bps to skip until the next)
    MM_m=split_mm[1] #contains 5mC mm tag, positions of methylations (how many bps to skip until the next)


    #add semicolon at the end
    MM_h=MM_h + ";"
    MM_m=MM_m + ";"
    #print(MM_h)
    #print(MM_m)

    ml_tags=read.get_tag("ML")
    length_of_array=int(len(ml_tags))
    half_length=int(length_of_array/2)

    ML_h=ml_tags[0:half_length] #contains probabilities of 5hmC methylation
    ML_m=ml_tags[half_length:length_of_array] #contains probabilities of 5mC methylation

    #print(ML_h)
    #print(ML_m)
    #print("\n")

    #NOW SET READ TAGS TO THESE NEW VALUES AND WRITE TO A FILE
    #5hmC
    read.set_tag("MM",MM_h)
    read.set_tag("ML",ML_h)
    output_5hmC.write(read)

    #5mC
    read.set_tag("MM",MM_m)
    read.set_tag("ML",ML_m)
    output_5mC.write(read)

samfile.close()
print("Done. Output files created.")
