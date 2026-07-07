CD45Expr_Comparison = function(genename){row_info <- as.data.frame(rowData(data_se))
target_row <- rownames(row_info[row_info$gene_name == genename, ])
gene_expr <- expression[target_row, ] #use expression referencing normalized data summary experiment
gene_name_group = paste(genename, "_group", sep = "")
clinical[[genename]] <- as.numeric(gene_expr[match(clinical$barcode, names(gene_expr))])
clinical<- clinical%>% arrange(clinical[[genename]])
clinical<- filter(clinical, shortLetterCode == "TP")
clinical[[gene_name_group]] <- case_when(clinical[[genename]] > quantile(clinical[[genename]], prob = 0.75, na.rm = TRUE) ~ "High",
                                         clinical[[genename]] < quantile(clinical[[genename]], prob = 0.25, na.rm = TRUE) ~ "Low",
) 


clean_df <- clinical[!is.na(clinical$time) & !is.na(clinical[[gene_name_group]]) & clinical$shortLetterCode %in% "TP", ]
target_row <- rownames(row_info[row_info$gene_name == "PTPRC" , ])
gene_expr2 <- expression[target_row, ]
clean_df$PTPRC <- as.numeric(gene_expr2[match(clean_df$barcode, names(gene_expr2))])
clean_df2 = select(clean_df, gene_name_group, genename, PTPRC)
clean_df2 = pivot_longer(clean_df2, c(genename, "PTPRC"), names_to = "Genes", values_to = "Expression")

wilcoxgroup <- paste("Expression ~ ", gene_name_group, sep = "")

stat.test <- compare_means(as.formula(wilcoxgroup), data = clean_df2, group.by = "Genes") %>% adjust_pvalue() %>% add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(data = clean_df2, formula = as.formula(wilcoxgroup), x = gene_name_group)
stat.test$p.adj.format <- p_format(stat.test$p.adj, accuracy = 0.01, leading.zero = FALSE)
bxplot <- ggboxplot(clean_df2, x = gene_name_group, y = "Expression", fill= "#E57373", facet.by="Genes") +  scale_y_log10(labels = scales::label_comma()) 

result = bxplot + stat_pvalue_manual(stat.test, label = "p.format", y.position = 5)
return(result)

}

CD45combined_boxplots <- map(GeneNames, CD45Expr_Comparison)


destination = "C:\\Users\\weniw\\Downloads\\CD45combined_boxplots.pdf"
pdf(file = destination)
par(mfrow = c(2,2))
CD45combined_boxplots
dev.off()