

GeneExpression = function(genename) { 
  row_info <- as.data.frame(rowData(data_se))
  target_row <- rownames(row_info[row_info$gene_name == genename, ])
  gene_expr <- expression[target_row, ]
  gene_name_expr = paste(genename, "_expr", sep = "")
  gene_name_group = paste(genename, "_group", sep = "")
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
  sf <- paste("Surv(time, status) ~", gene_name_group)
  fit <- surv_fit(as.formula(sf), data = clean_df)
  
  result = ggsurvplot(fit, 
                      data = clean_df,
                      pval = TRUE, 
                      risk.table = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      title = paste("Overall Survival", genename, "Expression"),
                      legend.labs = c(paste("High", genename), paste("Low", genename)))
  
  pvalue_result = surv_pvalue(fit)
  
  return_list = list(plot=result, p_value = pvalue_result)
  return(return_list)
}

