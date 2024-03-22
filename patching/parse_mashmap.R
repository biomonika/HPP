#!/usr/bin/env Rscript
#this script parses mashmap output in order to find out which breakpoints can be pursued for patching

library(dplyr)
options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No filename provided.")
}

file <- args[1]

mashmap<-as.data.frame(read.table(file,skip=0,sep=" "))
colnames(mashmap)<-c("contig","contig_length","cstart","cend","strand","chr","chr_length","chstart","chend","identity")

chr_value<-as.character(unique(sort(mashmap$chr)))

mashmap<-mashmap[
  with(mashmap, order(chstart)),
]
mashmap <- mashmap %>%
  filter(contig_length > 1000000)
print(mashmap)

# Function to check if coordinates are nested within any other row
is_nested <- function(index, df) {
  status_start<-df$chstart < df$chstart[index]
  status_end<-df$chend > df$chend[index]
  status<-(status_start & status_end)
  number_of_hits<-max(0,as.numeric(table(status)["TRUE"]),na.rm=TRUE)
  nested<-(number_of_hits>0)
  #print(paste(index,nested))
  return(nested)
}

# Function to filter rows based on nested coordinates
filter_nested_coordinates <- function(df) {
  filtered_rows <- sapply(1:nrow(df), function(i) !is_nested(i, df))
  filtered_df<-df[filtered_rows, ]
  filtered_df <- filtered_df %>% mutate(gap_size = lead(chstart) - chend)
  #filtered_df <- filtered_df %>% mutate(feasible = ifelse((gap_size <= 1000000) & (gap_size > 0) & !is.na(gap_size), TRUE, FALSE))
  return(filtered_df)
}

  filtered_chromosome_data <- filter_nested_coordinates(mashmap)
  print(filtered_chromosome_data)
  
  if (nrow(filtered_chromosome_data) > 1) {
    for (i in 1:(nrow(filtered_chromosome_data) - 1)) {
        print("Found neighboring contigs.")
        neighboring_contigs <- filtered_chromosome_data[i:(i + 1),]
        first <- as.data.frame(neighboring_contigs[1,])
        second <- as.data.frame(neighboring_contigs[2,])
        flanklength <- 1000000
        
        first_entry_bed <-
          paste(
            as.character(first$contig),
            as.numeric(first$cend) - flanklength,
            first$cend,
            as.character(chr_value),
            sep = '\t'
          )
        second_entry_bed <-
          paste(
            as.character(second$contig),
            0,
            as.numeric(0 + flanklength),
            as.character(chr_value),
            sep = '\t'
          )
        
        write.table(
          as.data.frame(rbind(
            first_entry_bed, second_entry_bed
          )),
          file = paste0(file,".order",i,".bed"),
          row.names = FALSE,
          col.names = FALSE,
          quote = FALSE
        )
    }
  } 
