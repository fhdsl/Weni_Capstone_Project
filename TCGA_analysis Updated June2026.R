if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
BiocManager::install("singscore")

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

#2. Extract Clinical and Expression data
clinical <- as.data.frame(colData(data_se))
clinical_flat <- as.data.frame(lapply(clinical, function(x) {
  if (is.list(x)) {
    return(sapply(x, paste, collapse = "; "))
  } else {
    return(x)
  }
}))
expression <- assay(data_se) # Rows = Genes, Cols = Patients; original code

#3. insert normalization step (new)
dds <- DESeqDataSet(data_se, design = ~1)
dds <- estimateSizeFactors(dds) 
normalized_counts <- counts(dds, normalized = TRUE)


expression <- normalized_counts #normalized_counts is the matrix 




row_info <- as.data.frame(rowData(data_se))
target_row <- rownames(row_info[row_info$gene_name == gene_name, ])
gene_expr <- expression[target_row, ]
gene_name_expr = paste(gene_name, "_expr", sep = "")
gene_name_group = paste(gene_name, "_group", sep = "")
clinical[[gene_name_expr]] <- as.numeric(gene_expr[match(clinical$barcode, names(gene_expr))])
clinical<- clinical%>% arrange(clinical[[gene_name_expr]])
clinical<- filter(clinical, shortLetterCode == "TP")

clinical[[gene_name_group]] <- case_when(clinical[[gene_name_expr]] > quantile(clinical[[gene_name_expr]], prob = 0.75, na.rm = TRUE) ~ "High",
                                          clinical[[gene_name_expr]] < quantile(clinical[[gene_name_expr]], prob = 0.25, na.rm = TRUE) ~ "Low",
) 
clinical$time <- ifelse(clinical$vital_status == "Dead", 
                         clinical$days_to_death, clinical$days_to_last_follow_up)
clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)

clean_df <- clinical[!is.na(clinical$time) & !is.na(clinical[[gene_name_group]]) & clinical$shortLetterCode %in% "TP", ]
f <- paste("Surv(time, status) ~", gene_name_group)
fit <- survfit(as.formula(f), data = clean_df)

ggsurvplot(fit, 
           data = clean_df,
           pval = TRUE, 
           risk.table = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           title = paste("Overall Survival", gene_name, "Expression"),
           legend.labs = c(paste("High", gene_name), paste("Low", gene_name)))

pvalue_result = surv_pvalue(fit)
