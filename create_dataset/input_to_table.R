library(R.utils)
read_counts_for_dir <- function(directory) {
  file_names<-list.files(path=directory, pattern='*.htseq.counts', recursive = TRUE)
  print(file_names)
  
  print(paste("    Reading file:", file_names[1]))
  counts_table <- read.table(paste(directory, file_names[1], sep=""), row.names=1)
  single_file_name <- unlist(strsplit(file_names[1], split='/', fixed=TRUE))[2]
  names(counts_table) <- c(single_file_name)
  
  for(i in 2:length(file_names)) {
    print(paste("    Reading file:", file_names[i]))
    file <- read.table(paste(directory, file_names[i], sep=""), row.names=1)
    
    print(unlist(strsplit(file_names[i], split='/', fixed=TRUE))[2])
    single_file_name <- unlist(strsplit(file_names[i], split='/', fixed=TRUE))[2]
    counts_table <- cbind(counts_table, file[1])
    names(counts_table)[names(counts_table) == "V2"] <- single_file_name
  }
  
  return(counts_table)
}

args <- commandArgs(trailingOnly = TRUE)
no_test_only_good_output <- read_counts_for_dir("/home/gebruiker/dataset/")
