library(edgeR)
library(R.utils)
library(ggplot2)

voom_normalize <- function(counts_table, design){
  #counts <- read.delim(counts_table)
  
  d0 <- DGEList(counts_table)
  d0 <- calcNormFactors(d0 , method =c("TMM"))
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,]
  voom(d, design = design, plot = T)
  v <- voom(d, design)
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

thyroid_counts_table <- read_counts_for_dir("/home/bram/Downloads/Thyroid/")

design_tsv <- read.table(file = '/home/bram/Documents/GenderBiasInHormonalCancers/design.csv', sep = '\t', header = TRUE)

sample_names <- colnames(thyroid_counts_table)
groups <- design_tsv[sample_names, ]
design <- model.matrix(~0+groups)

normalized_thyroid_counts <- voom_normalize(thyroid_counts_table, design)

fit <- lmFit(normalized_thyroid_counts, design)

contr <- makeContrasts(groupsF - groupsM, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))

gene <- "ENSG00000241859.5"
counts_for_top_gene <- as.data.frame(normalized_thyroid_counts$E[gene,])

rownames_design <- rownames(design_tsv)
female_htseq_files <- subset(design_tsv, Condition=="F")
male_htseq_files <- subset(design_tsv, Condition=="M")

counts_for_top_gene_f <- counts_for_top_gene[rownames(female_htseq_files), ]
counts_for_top_gene_m <- counts_for_top_gene[rownames(male_htseq_files), ]

counts_for_top_gene_f <- counts_for_top_gene_f[!is.na(counts_for_top_gene_f)]
counts_for_top_gene_m <- counts_for_top_gene_m[!is.na(counts_for_top_gene_m)]

df <- as.data.frame(matrix(nrow = 34, ncol = 2))
colnames(df) <- c("Male", "Female")

df$Female <- counts_for_top_gene_f
df$Male <- NA
df$Male[1:length(counts_for_top_gene_m)] <- counts_for_top_gene_m


boxplot(df, main = paste("Ensembl gene: ", gene, "\n adj.p.value: ", top.table[gene, ]$adj.P.Val, "\n logFC: ", top.table[gene, ]$logFC),
        ylab = "Normalized gene counts (voom)", xlab = "Gender")
