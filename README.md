# UCEC_RNASeq
# TCGA-UCEC RNA-Seq Analysis

## Project Description
This R-based workflow analyzes RNA-Seq data from TCGA-UCEC (Uterine Corpus Endometrial Carcinoma) project to identify differentially expressed genes between tumor and normal samples.

## Prerequisites
- R (version 4.0 or higher)
- RStudio
- Required R packages:
  - BiocManager
  - TCGAbiolinks
  - tidyverse
  - maftools
  - pheatmap
  - SummarizedExperiment

## Installation
```R
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
install.packages("tidyverse")
BiocManager::install("maftools")
install.packages("pheatmap")
BiocManager::install("SummarizedExperiment")
BiocManager::install("TCGAbiolinks")

## Workflow Steps

## Data Retrieval

## Query and download TCGA-UCEC RNA-Seq data
## Filter and prepare expression data


## Differential Expression Analysis

## Normalize gene expression data
## Identify differentially expressed genes
## Filter significant DEGs


## Results

Generate lists of upregulated and downregulated genes
Create visualization plots
Perform enrichment analysis



## File Structure

dataDEGs.csv: All differentially expressed genes
significant_degs.csv: Significantly differentially expressed genes
upregulated_genes.csv: Upregulated genes
downregulated_genes.csv: Downregulated genes
enrichment_barplot.pdf: GO enrichment visualization

## Usage

Set your working directory
Run the analysis scripts in sequence
Check output files in the results directory


TCGA Research Network
BioConductor Community
