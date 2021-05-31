#!/usr/bin/env Rscript
library(R.utils)
library(edgeR)
library(limma)
library(ggplot2)
library(readr)

# rawCountTable <- read.table('merged_files.csv', header = TRUE, sep = '\t', row.names = 1)
# 
# #rawCountTable <- rawCountTable[-c(1,2), -c(1,2,3,4,5)] 
# 
# rawCountTable <- as.matrix(rawCountTable)

normalization <- function(counts_table, design){
  # Make DGE list of the counts_label --> merged htseq counts files
  d0 <- DGEList(counts_table)
  # normalize the DGE list using the TMM method
  d0 <- calcNormFactors(d0 , method =c("TMM"))
  # do a cutoff
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,]
  # use voom to normalize 
  # voom(d, design = design, plot = T) --> used to make voom plots
  print(dim(d))
  print(dim(design))
  v <- voom(d, design)
  # return the voom data
  return(v)
}


# args[1] is the design.csv from the preprocessing and args[2] is file with the merged htseq counst file
args <- commandArgs(trailingOnly = TRUE)

# rewrite the design to a nummeric design
design_tsv <- read.table(file = args[1], sep = '\t', header = TRUE)

merged_files <- read_csv(args[2])
gene_row_names <- merged_files$X1
merged_files$X1 <- NULL
rownames(merged_files) <- gene_row_names

print(design_tsv)
sample_names <- colnames(merged_files)
print(sample_names)
groups <- design_tsv[sample_names, ]
print(groups)
design <- model.matrix(~0+groups)
print(design)

## design_tsv <- read.table(file = args[1], sep = '\t', header = TRUE)
# design <- read.table(file = args[1], sep = ',', header = TRUE)

# voom DGE data
voomDGE <- normalization(merged_files, design)

save(voomDGE, file = "voomDGE.RData")

# example terminal command 
# Rscript ./normalization.R '/home/gebruiker/design.csv' '/home/gebruiker/dataset/Cervix/'