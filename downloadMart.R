# This script is meant to download the emsembl dataset to the local machine,
# so that dataset updates do not make pathways dissappear
# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(biomaRt)

# Connect to the Ensembl biomart database
ensembl_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl")

mart_attributes <- listAttributes(ensembl_mart)

# Save the database to the local machine
save(ensembl_mart, file = "mmulatta_gene_ensembl.RData")
