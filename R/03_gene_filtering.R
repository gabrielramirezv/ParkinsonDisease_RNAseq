# Import libraries
library(SummarizedExperiment))

# Check quality of the data
rse_gene_SRP194918$assigned_gene_prop <-
  rse_gene_SRP194918$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP194918$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP194918$assigned_gene_prop)

# Observe differences among samples injected with different substances
with(colData(rse_gene_SRP194918),
     tapply(assigned_gene_prop,
            sra_attribute.injected_with,
            summary))

# Plot the distribution of the assigned gene proportion
hist(rse_gene_SRP194918$assigned_gene_prop, col = "purple")

# Save the unfiltered data
rse_gene_SRP194918_unfiltered <- rse_gene_SRP194918

# Filter genes with low expression
rse_gene_SRP194918 <-
  rse_gene_SRP194918[edgeR::filterByExpr(rse_gene_SRP194918),]

# Plot the distribution of the assigned gene proportion after filtering
hist(rse_gene_SRP194918$assigned_gene_prop, col = "purple")

# Check the expression of genes that were not filtered out
gene_means <- rowMeans(assay(rse_gene_SRP194918, "counts"))
summary(gene_means)

# Check the percentage of genes that were not filtered out
round(nrow(rse_gene_SRP194918) / nrow(rse_gene_SRP194918_unfiltered) * 100, 2)
