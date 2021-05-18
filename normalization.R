normalization <- function(counts_table, design){
  library(edgeR)
  #counts <- read.delim(counts_table)
  
  d0 <- DGEList(counts_table)
  d0 <- calcNormFactors(d0 , method =c("TMM"))
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,]
  d
  voom(d, design = design, plot = T)
  v <- voom(dge, design)
  return(v)
}
