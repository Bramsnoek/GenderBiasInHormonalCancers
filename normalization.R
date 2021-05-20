#!/usr/bin/env Rscript
library(R.utils)
library(edgeR)

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
  v <- voom(dge, design)
  # return the voom data
  return(v)
}


# args[1] is the design.csv from the preprocessing and args[2] is file with the merged htseq counst file
args <- commandArgs(trailingOnly = TRUE)

# rewrite the design to a nummeric design 
design_tsv <- read.table(file = args[1], sep = '\t', header = TRUE)
sample_names <- colnames(args[2])
groups <- design_tsv[sample_names, ]
design <- model.matrix(~0+groups)

# voom DGE data
voomDGE <- normalization(args[1], design)

save(voomDGE, file = "voomDGE.RData")

# example terminal command 
# Rscript ./normalization.R '/home/gebruiker/design.csv' '/home/gebruiker/dataset/Cervix/'