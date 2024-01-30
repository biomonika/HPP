library("gplots")

#USE THE APPROPRIATE INPUT FILE 
input_file<-"summary.pairwise.minimap2.PAN027.diploid.complete.poreC.complete.txt"
data <- read.table(input_file, header = FALSE, col.names = c("Value", "Column1", "Column2"))

#MODIFY ROWS AND COLUMNS TO MATCH YOUR DATA
#p-arms
col_names_PAN027<-c("haplotype1-0000002","haplotype2-0000072","unassigned-0001023","unassigned-0001559","unassigned-0001128","unassigned-0001092","unassigned-0001620","unassigned-0000644")
#q-arms
row_names_PAN027<-c("haplotype1-0000008","haplotype1-0000009","haplotype1-0000015","haplotype1-0000019","haplotype2-0000070","haplotype2-0000071","haplotype2-0000087","unassigned-0000703")

col_names<-col_names_PAN027
row_names<-row_names_PAN027

# Create an empty data frame
heatmap_df <- data.frame(matrix(NA, nrow = length(row_names), ncol = length(col_names)))
rownames(heatmap_df) <- row_names
colnames(heatmap_df) <- col_names

# Fill in the data frame with the values
for (i in 1:length(col_names)) {
  for (j in 1:length(row_names)) {
  chromosome<-col_names[i]
  unassigned<-row_names[j]
  #print(paste(chromosome,unassigned))
  entry1<-data[data$Column1==chromosome & data$Column2==unassigned,]
  entry2<-data[data$Column1==unassigned & data$Column2==chromosome,]
  
  contacts1_chromosome<-data[data$Column1==chromosome,] #all contacts for a chromosome from Column1
  contacts2_chromosome<-data[data$Column2==chromosome,] #all contacts for a chromosome from Column2
  sum_chromosomal_contacts<-sum(contacts1_chromosome$Value)+sum(contacts2_chromosome$Value)
  mean_chromosomal_contacts<-mean(c(contacts1_chromosome$Value,contacts2_chromosome$Value))
  
  contacts1_unassigned<-data[data$Column1==unassigned,] #all contacts for an unassigned from Column1
  contacts2_unassigned<-data[data$Column2==unassigned,] #all contacts for an unassigned from Column2
  sum_unassigned_contacts<-sum(contacts1_unassigned$Value)+sum(contacts2_unassigned$Value)
  
  entry1_value <- if (nrow(entry1) == 0) 0 else sum(entry1$Value, na.rm = TRUE)
  entry2_value <- if (nrow(entry2) == 0) 0 else sum(entry2$Value, na.rm = TRUE)
  
  entry<-entry1_value+entry2_value
  
  proportion_chromosome<-round(entry/sum_chromosomal_contacts*100)
  proportion_unassigned<-round(entry/sum_unassigned_contacts*100)
  mean_unassigned<-round(entry/mean_chromosomal_contacts*100)
  

  heatmap_df[unassigned, chromosome] <- ifelse(is.na(entry), 0, mean_unassigned)
  }
}

heatmap.2(
  as.matrix(heatmap_df),
  Rowv = FALSE,
  Colv = FALSE,
  dendrogram = "none",
  cellnote = heatmap_df,
  notecol = "black",
  col = rev(heat.colors(12)),
  main = "Heatmap",
  xlab = "",
  ylab = "",
  key = FALSE,
  key.title = "Value",
  trace = "none",
  cexRow = 0.6, 
  cexCol = 0.6,
  keysize=1,
)
