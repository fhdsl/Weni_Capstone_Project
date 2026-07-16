#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.22")
#BiocManager::install("TCGAbiolinks")
#BiocManager::install("DESeq2")
#BiocManager::install("singscore")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(singscore)
library(purrr)
library(ggpubr)
library(rstatix)
library(paletteer)
library(scales)


#1. Query for RNA-Seq data
query <- GDCquery(
  project = "TCGA-KIRC", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data_se <- GDCprepare(query)


#3. Normalize data
dds <- DESeqDataSet(data_se, design = ~1)
dds <- estimateSizeFactors(dds) 
normalized_counts <- counts(dds, normalized = TRUE)

#Analysis 1: Survival plots for a select list of genes.
source("TCGA_analysis Function 2026.R")
GeneNames = c("ADAMTS2", "ADAMTSL4", "ANGPTL4", "ASPN", "BMP1", "CD109", "COL12A1", "COL14A1", "COL16A1", "COL1A1", "COL1A2", "COL5A2", "COL6A1", "COL6A2", "COL8A1", "CRLF1", "DPT", "EMILIN1", "FBN1", "FN1", "LOX", "LOXL2", "LTBP1", "MFAP2", "MFAP5", "MMP19", "POSTN", "SERPINE1", "TGFBI", "THBS1", "THBS2", "THSD4", "TIMP1", "TNC","TSKU","VWA1")
results <- map(GeneNames, SingleGeneSurvivalAnalysis, data_se, normalized_counts)
pdf(file = "SingleGeneSurvivalAnalysis.pdf")
map(results, 1)
dev.off()

#Analysis 2: CD8 boxplots
source("CD8_TCGA Function.R")
combined_boxplots <- map(GeneNames, CD8Expr_Comparison, data_se, normalized_counts)
pdf(file = "combined_boxplots.pdf")
par(mfrow = c(2,2))
combined_boxplots
dev.off()

#Analysis 3: CD45 boxplots
source("CD45_TCGA Function.R")
CD45combined_boxplots <- map(GeneNames, CD45Expr_Comparison, data_se, normalized_counts)
pdf(file = "CD45combined_boxplots.pdf")
par(mfrow = c(2,2))
CD45combined_boxplots
dev.off()

#Analysis 4: CD68 boxplots
source("CD68_TCGA Function.R")
CD68combined_boxplots <- map(GeneNames, MacroExpr_Comparison, data_se, normalized_counts)
pdf(file = "CD68combined_boxplots.pdf")
par(mfrow = c(2,2))
CD68combined_boxplots
dev.off()


#Analysis 5: Composite Gene Score Survival Analysis 
source("Updated Composite Score and OS.R")
composite_analysis <- compositeSurvivalAnalysis(data_se, normalized_counts)
pdf(file = "CompositeSurvivalAnalysis.pdf")
composite_analysis$plot
dev.off()
