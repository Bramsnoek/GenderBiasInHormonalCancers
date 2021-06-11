#!/usr/bin/env Rscript
library(R.utils)
library(edgeR)
library(limma)
library(ggplot2)
library(readr)

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
  v <- voom(d, design)
  # return the voom data
  return(v)
}

linear_model_de <- function(normalized_counts, design) {
  fit <- lmFit(normalized_counts, design)
  contr <- makeContrasts(groupsF - groupsM, levels = colnames(coef(fit)))
  
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  qvals <- p.adjust(tmp$p.value, method = 'BH')
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
}

# args[1] is the design.csv from the preprocessing and args[2] is file with the merged htseq count file, args[3] is the output dir
args <- commandArgs(trailingOnly = TRUE)

# rewrite the design to a numeric design
design_tsv <- read.table(file = args[1], sep = '\t', header = TRUE)

merged_files <- read_csv(args[2])
gene_row_names <- merged_files$X1
merged_files$X1 <- NULL
rownames(merged_files) <- gene_row_names

sample_names <- colnames(merged_files)
groups <- design_tsv[sample_names, ]
design <- model.matrix(~0+groups)

# voom DGE data
voomDGE <- normalization(merged_files, design)

lm_results <- linear_model_de(voomDGE, design)
lm_results_0_0_5 <- lm_results[lm_results$adj.P.Val < 0.05, ]

write.csv(lm_results, paste(args[3], 'csv/lm_results.csv'), row.names = TRUE)
write.csv(lm_results_0_0_5, paste(args[3], 'csv/lm_results_005.csv'), row.names = TRUE)

# example terminal command 
# Rscript ./normalization.R '/home/gebruiker/design.csv' '/home/gebruiker/dataset/Cervix/' '/home/gebruiker/de_output/'
