#!/usr/bin/env Rscript
library(R.utils)
read_counts_for_dir <- function(directory) {
  file_names<-list.files(path=directory, pattern='*.htseq.counts')
  print(file_names)
  
  print(paste("    Reading file:", file_names[1]))
  counts_table <- read.table(paste(directory, file_names[1], sep=""), row.names=1)
  names(counts_table) <- c(paste(directory, file_names[1], sep=""))
  
  for(i in 2:length(file_names)) {
    print(paste("    Reading file:", file_names[i]))
    file <- read.table(paste(directory, file_names[i], sep=""), row.names=1)
    counts_table <- cbind(counts_table, file[1])
    names(counts_table)[names(counts_table) == "V2"] <- paste(directory, file_names[i], sep="")
  }
  
  return(counts_table)
}
# args[1] is the path to the input.tsv and args[2] is the path to the map with the hseqcounts which you want to use as input
args <- commandArgs(trailingOnly = TRUE)

# Making the design on the basis of the input.tsv
design <- read.table(file = args[1], sep = '\t', header = TRUE)

design[,1]
sort(row.names(design[,1]))
rownames(design) <-paste(args[2], design[,1], sep="")  
design[,1] <- NULL

# all files merged
merged_files <- read_counts_for_dir(args[2])

# example terminal command 
# Rscript ./input_to_table.R '/home/gebruiker/input.tsv' '/home/gebruiker/dataset/Cervix/'
