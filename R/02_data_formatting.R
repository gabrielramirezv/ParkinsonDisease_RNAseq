# Import libraries
library(recount3)
library(SummarizedExperiment)

# Transfrom nucleotide counts to read counts
assay(rse_gene_SRP194918, "counts") <- compute_read_counts(rse_gene_SRP194918)

# Add important data for analysis
rse_gene_SRP194918 <- expand_sra_attributes(rse_gene_SRP194918)

# Check the data format
rse_gene_SRP194918$sra.sample_attributes[1:10]

# Review data of interest
colData(rse_gene_SRP194918)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP194918)))
]

# Convert data types to factors
rse_gene_SRP194918$sra_attribute.age <-
  factor(rse_gene_SRP194918$sra_attribute.age)
rse_gene_SRP194918$sra_attribute.injected_with <-
  factor(rse_gene_SRP194918$sra_attribute.injected_with)
rse_gene_SRP194918$sra_attribute.source_name <-
  factor(rse_gene_SRP194918$sra_attribute.source_name)
rse_gene_SRP194918$sra_attribute.strain <-
  factor(rse_gene_SRP194918$sra_attribute.strain)
rse_gene_SRP194918$sra_attribute.tissue <-
  factor(rse_gene_SRP194918$sra_attribute.tissue)
rse_gene_SRP194918$sra_attribute.treatment <-
  factor(rse_gene_SRP194918$sra_attribute.treatment)

# Summarize data
summary(as.data.frame(colData(rse_gene_SRP194918)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP194918)))
]))

## There are 17 samples
## While the only relevant difference is the injection, we will analyze the
## effects that it has in gene expression
