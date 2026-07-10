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