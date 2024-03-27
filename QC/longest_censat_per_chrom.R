#!/usr/bin/Rscript
# This script reads censat annotation track in bed format 
# and outputs the contig names and the names of the most  
# abundant satellite feature on that contig
# 
# Written by Monika Cechova (mcechova@ucsc.edu)

# Read the input filename from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Usage: Rscript longest_censat_per_chrom.R input_censat_filename")
}

# Read the input file
input_censat_file <- args[1]
censat<-read.table(input_censat_file, header = FALSE, skip = 1)
colnames(censat)<-c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb")

require(dplyr)

censat_by_chrom_and_name <- censat %>% group_by(chrom,name)
censat_by_chrom_and_name <- summarise(censat_by_chrom_and_name, total_length = sum(chromEnd - chromStart))
longest_censat_per_chrom <- censat_by_chrom_and_name %>% slice_max(total_length) %>%
  ungroup()

output_file <- gsub("\\.bed$", "_longest_censat.txt", input_censat_file)
write.table(longest_censat_per_chrom, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
