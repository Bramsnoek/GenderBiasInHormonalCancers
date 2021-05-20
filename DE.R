library(edgeR)
library(R.utils)

voom_normalize <- function(counts_table, design){
  #counts <- read.delim(counts_table)
  
  d0 <- DGEList(counts_table)
  d0 <- calcNormFactors(d0 , method =c("TMM"))
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,
  voom(d, design = design, plot = T)
  v <- voom(dge, design)
  return(v)
}

read_counts_for_dir <- function(directory) {
  file_names<-list.files(path=directory, pattern='*.htseq.counts')
  print(file_names)
  
  print(paste("    Reading file:", file_names[1]))
  counts_table <- read.table(paste(directory, file_names[1], sep=""), row.names=1)
  names(counts_table) <- c(paste(file_names[1], sep=""))
  
  for(i in 2:length(file_names)) {
    print(paste("    Reading file:", file_names[i]))
    file <- read.table(paste(directory, file_names[i], sep=""), row.names=1)
    counts_table <- cbind(counts_table, file[1])
    names(counts_table)[names(counts_table) == "V2"] <- paste(file_names[i], sep="")
  }
  
  return(counts_table)
}

thyroid_counts_table <- read_counts_for_dir('C:/Users/Bram Snoek/thyroid/')

design_tsv <- read.table(file = 'C:/Users/Bram Snoek/GenderBiasInHormonalCancersGit/design.csv', sep = '\t', header = TRUE)

sample_names <- colnames(thyroid_counts_table)
groups <- design_tsv[sample_names, ]
design <- model.matrix(~0+groups)

normalized_thyroid_counts <- voom_normalize(thyroid_counts_table, 'C:/Users/Bram Snoek/GenderBiasInHormonalCancersGit/design.csv')

