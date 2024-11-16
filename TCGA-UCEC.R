##### Section 1: Package Installation and Loading####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the latest version of BiocManager
install.packages("BiocManager")

# Load BiocManager
library(BiocManager)

# Install other required packages
install.packages("tidyverse")
BiocManager::install("maftools")
install.packages("pheatmap")
BiocManager::install("SummarizedExperiment")

# Load necessary libraries
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)

##### Section 2: Retrieve and Explore TCGA Data####

#Initial Data Exploration
gdcprojects <- getGDCprojects()
gdcprojects
getProjectSummary("TCGA-UCEC")

# Building a query for transcriptome profiling data
query_TCGA <- GDCquery(project = 'TCGA-UCEC',
                       data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)

# Section 3: Download and Prepare Gene Expression Data
query_TCGA <- GDCquery(project = 'TCGA-UCEC',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')

res <- getResults(query_TCGA)

summary(factor(res$sample_type))

# Set the destination directory
# Replace this path with your own directory path
destination <- "C:/Users/YOUR_USERNAME/TCGA_Data"  # or wherever you want to save

# Create directory if it doesn't exist
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

# Download data
GDCdownload(query_TCGA, dir = destination)


##### Section 3: Prepare the downloaded data####
# Step 1: Prepare RNA-seq data
# GDCprepare converts downloaded data into a SummarizedExperiment object
# This object contains: gene expression counts, sample metadata, and feature data
UCEC_Rnaseq_SE <- GDCprepare(query_TCGA)

# Step 2: Extract expression matrix
# assay() gets the raw count matrix from the SE object
# 'unstranded' is used because the RNA-seq protocol didn't preserve strand information
# Matrix dimensions: rows = genes, columns = samples
UCECMatrix <- assay(UCEC_Rnaseq_SE, "unstranded")

# Step 3: Quality control and preprocessing
# - Creates visualization plots for quality assessment
# - Identifies and flags potential outlier samples
# - Performs initial data cleaning steps
# - Returns preprocessed expression data for downstream analysis
UCEC.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(UCEC_Rnaseq_SE)


##### Section 4: Normalization of genes####
dataNorm <- TCGAanalyze_Normalization(
  tabDF = UCEC.RNAseq_CorOutliers, 
  geneInfo =  geneInfoHT
)

normalized_count_matrix <- dataNorm

# Write normalized_count_matrix to a CSV file in the specified directory
write.csv(normalized_count_matrix, 
          file = "C:/Users/YOUR_USERNAME/YOUR_FOLDER/normalized_counts.csv", 
          row.names = FALSE)


##### Section 5: Data Filtering####

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)
##### Section 6: Sample Selection####
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

# Print number of samples
print(paste("Number of normal samples:", length(samplesNT)))
print(paste("Number of tumor samples:", length(samplesTP)))

# Check first few samples of each type
print("First few normal samples:")
head(samplesNT)
print("First few tumor samples:")
head(samplesTP)


##### Section 7: Differential Expression Analysis####
# Diff.expr.analysis (DEA) 
#Compares tumor vs normal samples
#Uses FDR < 0.01 and logFC > 1 cutoffs

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

# Save the results
write.csv(dataDEGs, "dataDEGs.csv", row.names = TRUE)

# Then read it (if needed)
data <- read.csv("dataDEGs.csv")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)


# DEGs table with expression values
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)

# Set output directory
output_dir <- "C:/Users/YOUR_USERNAME/Desktop/Senior thesis"
dir.create(output_dir, showWarnings = FALSE)

# Save initial DEGs
write.csv(dataDEGs, file.path(output_dir, "dataDEGs.csv"), row.names = TRUE)

# Print summary statistics
num_DEGs <- nrow(dataDEGs)
print(paste("Total number of DEGs:", num_DEGs))

# Read and filter DEGs
degs <- read.csv(file.path(output_dir, "dataDEGs.csv"))
significant_degs <- subset(degs, logFC > 2 | logFC < -2)
print(paste("Number of significant DEGs:", nrow(significant_degs)))

# Separate up/down regulated genes
upregulated_genes <- subset(significant_degs, logFC > 2)
downregulated_genes <- subset(significant_degs, logFC < -2)

# Save results
write.csv(significant_degs, file.path(output_dir, "significant_degs.csv"), row.names = TRUE)
write.csv(upregulated_genes, file.path(output_dir, "upregulated_genes.csv"), row.names = TRUE)
write.csv(downregulated_genes, file.path(output_dir, "downregulated_genes.csv"), row.names = TRUE)

# Print summary
print(paste("Upregulated genes:", nrow(upregulated_genes)))
print(paste("Downregulated genes:", nrow(downregulated_genes)))

# Optional: View results
View(significant_degs)
###### Section 6: TCGAanalyze_EAcomplete & TCGAvisualize_EAbarplot: Enrichment Analysis#####

library(TCGAbiolinks)
# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP), 
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = Genelist, 
  nBar = 10
)




