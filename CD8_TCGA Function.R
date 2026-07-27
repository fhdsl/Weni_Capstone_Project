CD8Expr_Comparison = function(genename, summarized_experiment, normalized_expression){
  clinical <- as.data.frame(colData(summarized_experiment))
  row_info <- as.data.frame(rowData(summarized_experiment))
  target_row <- rownames(row_info[row_info$gene_name == genename, ])
  gene_expr <- normalized_expression[target_row, ] #use expression referencing normalized data summary experiment
  gene_name_group = paste(genename, "_group", sep = "")
  clinical[[genename]] <- as.numeric(gene_expr[match(clinical$barcode, names(gene_expr))])
  clinical<- clinical%>% arrange(clinical[[genename]])
  clinical<- filter(clinical, shortLetterCode == "TP")
  clinical[[gene_name_group]] <- case_when(clinical[[genename]] > quantile(clinical[[genename]], prob = 0.75, na.rm = TRUE) ~ "High",
                                           clinical[[genename]] < quantile(clinical[[genename]], prob = 0.25, na.rm = TRUE) ~ "Low") 
  clinical$time <- ifelse(clinical$vital_status == "Dead", 
                          clinical$days_to_death, clinical$days_to_last_follow_up)
  clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)
  
  clean_df <- clinical[!is.na(clinical[[gene_name_group]]) & clinical$shortLetterCode %in% "TP", ]
  target_row <- rownames(row_info[row_info$gene_name == "CD8A" , ])
  gene_expr2 <- normalized_expression[target_row, ]
  clean_df$CD8A <- as.numeric(gene_expr2[match(clean_df$barcode, names(gene_expr2))])
  target_row <- rownames(row_info[row_info$gene_name == "CD8B" , ])
  gene_expr3 <- normalized_expression[target_row, ]
  clean_df$CD8B <- as.numeric(gene_expr3[match(clean_df$barcode, names(gene_expr3))])
  
  clean_df2 = dplyr::select(clean_df, gene_name_group, genename, CD8A, CD8B)
  clean_df2 = pivot_longer(clean_df2, c(genename, "CD8A", "CD8B"), names_to = "Genes", values_to = "Expression")
  
  wilcoxgroup <- paste("Expression ~ ", gene_name_group, sep = "")
  
  stat.test <- compare_means(as.formula(wilcoxgroup), data = clean_df2, group.by = "Genes") %>% adjust_pvalue() %>% add_significance("p.adj")
  stat.test <- stat.test %>% add_xy_position(data = clean_df2, formula = as.formula(wilcoxgroup), x = gene_name_group)
  stat.test$p.adj.format <- p_format(stat.test$p.adj, accuracy = 0.01, leading.zero = FALSE)
  bxplot <- ggboxplot(clean_df2, x = gene_name_group, y = "Expression", fill= "#E57373", facet.by="Genes") +  scale_y_log10() 
  
  
  result = bxplot + stat_pvalue_manual(stat.test, label = "p.format", y.position = 4)
  return(result)

}
