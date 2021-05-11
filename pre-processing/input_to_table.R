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
args <- commandArgs(trailingOnly = TRUE)

input_tsv <- read.table(file = args[1], sep = '\t', header = TRUE)

input_tsv[,1]
sort(row.names(input_tsv[,1]))
rownames(input_tsv) <-paste(args[2], input_tsv[,1], sep="")  
input_tsv[,1] <- NULL

merged_files <- read_counts_for_dir(args[2])

