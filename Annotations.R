#Install and load required packages
if(!requireNamespace("BiocManager",quietly = TRUE))
install.packages("BiocManager")  

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "pathview"
))

#Load them
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
library(dplyr)

#use significant DEGs 
sig_final <- res_unique %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2)
genes_final <- sig_final$hgnc_symbol
genes_final <- unique(genes_final)
genes_final <- genes_final[genes_final != ""]
length(genes_final)

#convert gene symbols
gene_df_final <- bitr(
  genes_final,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
nrow(gene_df_final)

#Functional annotation
###biological process
ego_bp_final <- enrichGO(
  gene = gene_df_final$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
length(ego_bp_final@gene)

#molecular function
ego_mf_final <- enrichGO(
  gene          = gene_df_final$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

#Cellular component
ego_cc_final <- enrichGO(
  gene = gene_df_final$ENTREZID,
  OrgDb =  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

#Visualization
dotplot(ego_bp_final,
        showCategory = 15,
        font.size = 12
)+
  ggtitle("GO Biological Process Enrichment (Significant DEGs)") +
  theme_minimal()

dotplot(
  ego_cc_final,
  showCategory = 15,
  font.size = 12
) +
  ggtitle("GO Cellular Component Enrichment") +
  theme_minimal()


#KEGG Pathway enrichment
ekegg <- enrichKEGG(
  gene = gene_df_final$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(ekegg)
ekegg <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)
head(gene_df_final$ENTREZID)

#Visualzie kEGG pathways
dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment - Lung Cancer")

head(ekegg@result$ID)
pathview(
  gene.data = sig$logFC,
  pathway.id = ekegg@result$ID[1],
  species = "hsa"
)
#save results
write.csv(as.data.frame(ego_bp), "GO_BP_results.csv")
write.csv(as.data.frame(ego_mf), "GO_MF_results.csv")
write.csv(as.data.frame(ego_cc), "GO_CC_results.csv")
write.csv(as.data.frame(ekegg),  "KEGG_results.csv")
length(unique(sig_final$hgnc_symbol))   # ≈ 160
nrow(gene_df)                           # slightly less (normal)
ego_bp <- enrichGO(...)
length(ego_bp@gene)

setwd("C:/Users/fathi/Downloads/lung_cancer_project")
deg <- read.csv("deg_outputs/DEG_significant.csv")
top <- deg[deg$adj.P.Val < 0.01 & abs(deg$logFC) > 2, ]
genes <- unique(top$hgnc_symbol)
length(genes)
writeLines(genes, "deg_outputs/hub_genes/ppi_input_genes.txt")
