ecm_signature <- c("ADAMTS2", "ADAMTSL4", "ANGPTL4", "ASPN", "BMP1", "CD109",
                   "COL12A1", "COL14A1", "COL16A1", "COL1A1", "COL1A2", "COL5A2",
                   "COL6A1", "COL6A2", "COL8A1", "CRLF1", "DPT", "EMILIN1",
                   "FBN1", "FN1", "LOX", "LOXL2", "LTBP1", "MFAP2", "MFAP5",
                   "MMP19", "POSTN", "SERPINE1", "TGFBI", "THBS1", "THBS2",
                   "THSD4", "TIMP1", "TNC", "TSKU", "VWA1")

length(ecm_signature)  # should be 36

# check for missing genes before scoring
missing_genes <- setdiff(ecm_signature, rownames(expression))
missing_genes
sum(ecm_signature %in% rownames(expression))


rownames(expression) <- rowData(data_se)$gene_name
sum(duplicated(rownames(expression)))

rankData <- rankGenes(expression)
sum(ecm_signature %in% rownames(rankData))  # should now show close to 36
scoredResults <- simpleScore(rankData, upSet = ecm_signature)

#subset data from clinical df
clinicalsing <- clinical|> dplyr::select(days_to_last_follow_up, days_to_death, vital_status, shortLetterCode)

#combine dataframes (clinicalsing + scoredResults); need to create ensembl_id into column given 

clinicalsing$ensembl_id <- rownames(clinicalsing)
scoredResults$ensembl_id <- rownames(scoredResults)

clinicalsing1 <- full_join(scoredResults, clinicalsing, by = "ensembl_id")

#reorganizing column order

clinicalsing1 <- clinicalsing1 %>% relocate(ensembl_id, TotalScore, .before = TotalDispersion)

#creating OS curves of all composite scores for 36 genes across 614 samples 

clinicalsing1<- filter(clinicalsing1, shortLetterCode == "TP")
clinicalsing1$scoregroup <- case_when(clinicalsing1$TotalScore > quantile(clinicalsing1$TotalScore, prob = 0.75, na.rm = TRUE) ~ "High",
                                         clinicalsing1$TotalScore < quantile(clinicalsing1$TotalScore, prob = 0.25, na.rm = TRUE) ~ "Low",
) 
clinicalsing1$time <- ifelse(clinicalsing1$vital_status == "Dead", 
                        clinicalsing1$days_to_death, clinicalsing1$days_to_last_follow_up)
clinicalsing1$status <- ifelse(clinicalsing1$vital_status == "Dead", 1, 0)

clean_clinicalsing1 <- clinicalsing1[!is.na(clinicalsing1$time) & !is.na(clinicalsing1$scoregroup) & clinicalsing1$shortLetterCode %in% "TP", ]
fit <- survfit(as.formula("Surv(time, status) ~ scoregroup"), data = clean_clinicalsing1)

ggsurvplot(fit, 
           data = clinicalsing1,
           pval = TRUE, 
           risk.table = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           title = paste("Overall Survival for All 36 Genes Composite Score"),
           legend.labs = c(paste("High Score Group"), paste("Low Score Group")))

pvalue_result = surv_pvalue(fit)

