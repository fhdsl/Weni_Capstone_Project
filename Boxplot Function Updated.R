gene_name_expr_High = paste(gene_name, "_expr_High", sep = "")
gene_name_expr_Low = paste(gene_name, "_expr_Low", sep = "")


clean_df[[gene_name_expr_High]] <- case_when(clean_df[[gene_name_group]] == "High"~ clean_df[[gene_name_expr]])
clean_df[[gene_name_expr_Low]] <- case_when(clean_df[[gene_name_group]]== "Low"~ clean_df[[gene_name_expr]])
clean_df$CD8A_exprHigh <- case_when(clean_df[[gene_name_group]] == "High"~ clean_df$CD8A_expr)
clean_df$CD8A_exprLow <- case_when(clean_df[[gene_name_group]] == "Low" ~ clean_df$CD8A_expr)
clean_df$CD8B_exprHigh <- case_when(clean_df[[gene_name_group]] == "High" ~ clean_df$CD8B_expr)
clean_df$CD8B_exprLow <- case_when(clean_df[[gene_name_group]] == "Low" ~ clean_df$CD8B_expr)





plotselect_df <- select(clean_df, [[gene_name_expr_High]], [[gene_name_expr_Low]], CD8A_expr_High, CD8A_expr_Low, CD8B_expr_High, CD8B_expr_Low)
boxplot_df <- pivot_longer(plotselect_df, everything())
cleanboxplot_df <- boxplot_df %>% drop_na()



cleanboxplot_df$name <- factor(cleanboxplot_df$name, levels=c("[[gene_name_expr]]_High", "CD8A_expr_High", "CD8B_expr_High", "[[gene_name_expr]]_Low", "CD8A_expr_Low", "CD8B_expr_Low"))


ggplot(cleanboxplot_df, aes(x=name, y=value, fill=name)) + 
  geom_boxplot(alpha=0.3) +
  scale_y_log10() +
  theme(legend.position="none") 