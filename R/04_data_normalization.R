# Import libraries
library(edgeR)

# Create the list for normalization
dge <- DGEList(counts = assay(rse_gene_SRP194918, "counts"),
               genes = rowData(rse_gene_SRP194918))

# Calculate the normalization factors
dge <- calcNormFactors(dge)
