#!/usr/bin/env Rscript
#This script classifies the assembly into maternal and paternal contigs using parental k-mers from PAN010 and PAN011
#Written by Monika Cechova mcechova@ucsc.edu

# Check if exactly two arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Usage: script.R maternal_file paternal_file")
}

# Extract filenames from command-line arguments
maternal_file <- commandArgs(trailingOnly = TRUE)[1]
paternal_file <- commandArgs(trailingOnly = TRUE)[2]

# Read maternal data
maternal <- read.table(maternal_file)
colnames(maternal) <- c("haplotype", "maternal")

# Read paternal data
paternal <- read.table(paternal_file)
colnames(paternal) <- c("haplotype", "paternal")

data<-merge(maternal,paternal)

#Normalize per coverage
data$maternal_normalized<-data$maternal/sum(data$maternal)*100
data$paternal_normalized<-data$paternal/sum(data$paternal)*100

#Calculate the ratio
data$ratio<-data$maternal_normalized-data$paternal_normalized

# Perform clustering
kmeans_result <- kmeans(data$ratio, centers = 2)

# Assign labels based on cluster centroids
centroid_values <- kmeans_result$centers
maternal_label <- which.max(centroid_values)
paternal_label <- which.min(centroid_values)

# Assign labels to haplotypes based on clustering results
data$parent_of_origin <- ifelse(kmeans_result$cluster == maternal_label, "maternal", "paternal")
data<-data[,c("haplotype","parent_of_origin")]

output_name<-gsub("maternal", "parental", maternal_file)
output_name<-gsub(".txt", "", output_name)

# Write the resulting table to a file
write.table(data, paste0("parent_of_origin.",output_name,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
