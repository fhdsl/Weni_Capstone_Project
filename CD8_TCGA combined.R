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
target_row <- rownames(row_info[row_info$gene_name == "CD8A" , ])
gene_expr2 <- expression[target_row, ]
clean_df$CD8A_expr <- as.numeric(gene_expr2[match(clean_df$barcode, names(gene_expr2))])
target_row <- rownames(row_info[row_info$gene_name == "CD8B" , ])
gene_expr3 <- expression[target_row, ]
clean_df$CD8B_expr <- as.numeric(gene_expr3[match(clean_df$barcode, names(gene_expr3))])
