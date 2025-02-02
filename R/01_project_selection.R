# Import libraries
library(recount3)
library(edgeR)

# List all available projects for mouse
mouse_projects <- available_projects(organism = "mouse")
mouse_projects

# Extract the project information
proj_info <- subset(
  mouse_projects,
  project == "SRP194918" & project_type == "data_sources"
)

# Create a recount object for the SRP058181 project
rse_gene_SRP194918 <- create_rse(proj_info)
rse_gene_SRP194918

# Analyze the project
rowData(rse_gene_SRP194918)
rowRanges(rse_gene_SRP194918)
