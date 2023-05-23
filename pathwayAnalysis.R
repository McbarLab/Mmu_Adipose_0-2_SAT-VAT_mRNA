# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if(!require("org.Mmu.eg.db", quietly = TRUE))
  BiocManager::install("org.Mmu.eg.db")
if(!require("EnhancedVolcano", quietly = TRUE))
  BiocManager::install('EnhancedVolcano')
if(!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(tidyverse)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(EnhancedVolcano)
library(AnnotationDbi)
library(biomaRt)

# Local database is used to ensure pathway does not change due to updates
emsembl_mart <- load("mmulatta_gene_ensembl.RData")

import_dataset <- function(filename){
  rawDGE <- read_xlsx(filename)
  # The following line is for Mark Berres dataset
  curatedDGE <- rawDGE[,c(1,2,3,6,7)] %>% na.omit()
  colnames(curatedDGE) <- c("Ensembl", "Symbol", "logFC", "pval", "FDR")
  
  # Add a column named direction, to show whether this gene is UP-regulated
  # or DOWN-regulated
  # Take into consideration of the FDR. Do not simply look at the sign.
  curatedDGE <- curatedDGE %>% 
    mutate(direction = case_when(FDR < fdrCutOff & logFC > logFoldChange ~ "UP",
                                 FDR < fdrCutOff & logFC < logFoldChange ~ "DOWN",
                                 FDR >= fdrCutOff ~ "NS")) 
  return(curatedDGE)
}

volcano_plot <- function(curatedDGE, title){
  volcanoPlot <- 
    curatedDGE %>% EnhancedVolcano(
      lab = NA,
      # lab = curatedDGE$Symbol,
      x = 'logFC',
      # IMPORTANT: Do NOT put logFDR here, it takes negative log by default
      y = 'FDR', 
      xlab = expression(paste("Log"[2],"(FC)")),
      ylab = expression(paste("-log"[10],"(FDR)")),
      pCutoff = fdrCutOff,
      FCcutoff = logFoldChange,
      cutoffLineType = 'twodash',
      # col = c("#999999", "#ed1c24", "#662d91"),
      # pointSize = 3.0,
      labSize = 4.0,
      boxedLabels = TRUE,
      colAlpha = 4/5,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      colConnectors = 'grey50',
      # Put a legend. Good habit :)
      legendLabels = c("Not Sig", "LogFC", "-logFDR", "LogFC and -logFDR"),
      legendPosition = 'right',
      legendLabSize = 7,
      legendIconSize = 1.7,
      title = paste(title,"Volcano",sep=" "),
      subtitle = '',
      caption = '',
      max.overlaps = 100
    )
  ggsave(path = "./Volcano_Plots",
         filename = paste(title,"Volcano.pdf",sep=" "), 
         volcanoPlot)
}

GSEA_plot <- function(curatedDGE, title){
  GSEA_genes <- curatedDGE$logFC
  # GSEA can take in both Ensembl or Symbol
  names(GSEA_genes) <- curatedDGE$Ensembl
  # ENTREZ ids for more accurate result
  GSEA_id <- bitr(names(GSEA_genes), fromType = "ENSEMBL",
                  toType = "ENTREZID", OrgDb = org.Mmu.eg.db)
  
  # Remove duplicated entries
  dedup_id <- GSEA_id[!duplicated(GSEA_id[c("ENSEMBL")]), ]
  dedup_id <- dedup_id[!duplicated(dedup_id[c("ENTREZID")]), ]
  
  # Remove genes that have more than 1 rows
  # Note: This removes ALL rows for those genes from the dataset
  # Extract these deduplicated entries from the original curatedDGE
  mappedIDs <- curatedDGE %>% 
    filter(!duplicated(Ensembl)) %>% 
    filter(Ensembl %in% dedup_id$ENSEMBL)
  mappedIDs$entrez <- dedup_id$ENTREZID
  
  # Further distill the data by removing nulls
  kegg_genes <- mappedIDs$logFC
  names(kegg_genes) <- mappedIDs$entrez
  kegg_genes <- kegg_genes %>% na.omit()
  kegg_genes <- sort(kegg_genes, decreasing = T)
  
  # Set the timeout limit for gseKEGG() and useMart()
  options(timeout = Inf)
  
  # This is the very step of running GSEA, takes VERY LONG TIME
  gsea_result <- gseKEGG(geneList = kegg_genes, 
                         organism = "mcc",
                         minGSSize = 4, 
                         maxGSSize = 500, 
                         pvalueCutoff = fdrCutOff, 
                         pAdjustMethod = "fdr",
                         keyType = "ncbi-geneid",
                         verbose = TRUE)
  
  # Output all pathways as text into a csv file
  # Translate the Entrez IDs into gene symbols for easier interpretation
  
  # Check if a dataframe is empty
  is_empty <- function(df) {
    return(nrow(df) == 0 && ncol(df) == 0)
  }
  # Get the gene symbols and add them to gsea_result@result
  if(is_empty(gsea_result)==TRUE){
    break
  }else{
    for (i in seq_along(gsea_result@result$core_enrichment)) {
      input_string = gsea_result@result$core_enrichment[i]
      input_ids = strsplit(input_string, "/")[[1]]
      output_ids = "external_gene_name"
      
      # Get the gene symbols
      gsea_genes = getBM(attributes = c(output_ids), 
                         filters = c("entrezgene_id"), 
                         values = input_ids, 
                         mart = ensembl_mart,
      )
      
      # Collapse the gene symbols into a single string separated by backslashes
      gene_string = paste(gsea_genes[,output_ids], collapse = "/")
      
      # Assign the gene symbol to the corresponding row of gsea_result@result
      gsea_result@result[i, "gene_symbol"] <- gene_string
    }
    
    # Write the dataframe into an external csv file
    write.csv(gsea_result@result,
              file = paste("GSEA_Pathway_csv/", title," gseaKEGG_top10.csv",sep=""),
              row.names = FALSE)
  }
  
  gsea_dotplot <- dotplot(gsea_result, 
                          showCategory = 10,
                          title = paste(title,"gseaKEGG_top10",sep=" "), 
                          split = ".sign") + facet_grid(.~.sign)
  
  ggsave(path = "./GSEA_Pathways",
         filename = paste(title,"gseaKEGG_top10.pdf",sep=" "), 
         gsea_dotplot)
}

# Write all DGE spreadsheet names into a list
DGE_folder_name <- "DGE"
DGE_list <- list.files(path = paste("./",DGE_folder_name,sep = ""))
DGE_count <- length(DGE_list)

logFoldChange = log2(1.5)
fdrCutOff = 0.05

# Iterate through all files to generate their figures in 1 click
for(DGE_index in 1:DGE_count){
  title <- sub(".xlsx*.","",
               sub(".*GENE_","",DGE_list[DGE_index]))
  curatedDGE <- import_dataset(paste(DGE_folder_name, DGE_list[DGE_index],sep="/"))
  volcano_plot(curatedDGE, title)
  GSEA_plot(curatedDGE, title)
}
