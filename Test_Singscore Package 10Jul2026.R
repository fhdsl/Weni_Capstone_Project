install.packages("devtools")
devtools::install_github('DavisLabooratory/singscore')
BiocManager::install("msigdb")
library(GSEABase)
tgfb_expr_10_se
colnames(tgfb_expr_10_se)
tgfb_gs_up
tgfb_gs_dn
length(GSEABase::geneIds(tgfb_gs_up))
length(GSEABase::geneIds(tgfb_gs_dn))
rankData <- rankGenes(tgfb_expr_10_se)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up, downSet = tgfb_gs_dn)

getStableGenes(5, type = 'carcinoma')

class(tgfb_gs_up)



up_genes = c("ADAMTS2", "ADAMTSL4", "ANGPTL4", "ASPN", "BMP1", "CD109", "COL12A1", "COL14A1", "COL16A1", "COL1A1", "COL1A2", "COL5A2", "COL6A1", "COL6A2", "COL8A1", "CRLF1", "DPT", "EMILIN1", "FBN1", "FN1", "LOX", "LOXL2", "LTBP1", "MFAP2", "MFAP5", "MMP19", "POSTN", "SERPINE1", "TGFBI", "THBS1", "THBS2", "THSD4", "TIMP1", "TNC","TSKU","VWA1")

custom_gene_set <- GeneSet(
  up_genes, 
  geneIdType = GenenameIdentifier(), # Specifies identifier type (e.g., Symbol, Entrez)
  setName = "TEST_Signature" # A unique name for your gene set
)



rankData = rankGenes(expression)

scoredf <- simpleScore(rankData, upSet = custom_gene_set, downSet = NULL)

#issue - singscore does not work because the expression matrix shows the genes in rows
#but gene names are in columns; hence it cannot find the gene names 

#scoredf <- simpleScore(rankData, upSet = up_genes, downSet = NULL)