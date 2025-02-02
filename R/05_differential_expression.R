# Import libraries
library(ggplot2)
library(limma)
library(pheatmap)
library(RColorBrewer)

# Plot gene expression according to the injection group
ggplot(as.data.frame(colData(rse_gene_SRP194918)),
       aes(y = assigned_gene_prop,
           x = sra_attribute.injected_with)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Injection Group") +
  theme_classic()

# Create the model matrix according to the injection group
mod <- model.matrix(~ sra_attribute.injected_with + assigned_gene_prop,
                    data = colData(rse_gene_SRP194918)
)
colnames(mod)
# Estimate the mean-variance relationship
vGene <- voom(dge, mod, plot = TRUE)

# Compute statistics by Bayes method
eb_results <- eBayes(lmFit(vGene))

# Extract the top genes
de_results <- topTable(
  eb_results,
  coef = c(2, 3),
  number = nrow(rse_gene_SRP194918),
  sort.by = "none"
)

# Calculate the number of significant genes
table(de_results$adj.P.Val < 0.05)

# Plot the results
plotMA(eb_results, coef = 2)
plotMA(eb_results, coef = 3)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
volcanoplot(eb_results, coef = 3, highlight = 3, names = de_results$gene_name)

# Review the top genes information
de_results[de_results$gene_name %in% c("Jchain", "Lox", "Slc2a5"), ]
de_results[de_results$gene_name %in% c("Gm15564", "Gm23935", "Mir6236"), ]

# Heatmap of the 50 top genes
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

# Create a dataframe with the injection group
df <- as.data.frame(colData(rse_gene_SRP194918)$sra_attribute.injected_with)

# Rename the columns and rows
colnames(df) <- c("Injection")
rownames(df) <- colnames(exprs_heatmap)

# Change the gene ID for the gene name
row.names(exprs_heatmap) <-
  de_results$gene_name[rank(de_results$adj.P.Val) <= 50]

# Create the heatmap
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row",
)

# MDS plot according to the injection group
col.group <- df$Injection
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(vGene$E, labels = df$Injection, col = col.group)
