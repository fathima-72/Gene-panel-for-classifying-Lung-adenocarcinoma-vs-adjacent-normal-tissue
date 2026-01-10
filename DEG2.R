# GSE40419 – limma DGE + biomaRt annotation
#############################################

library(limma)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(dplyr)

# ---------- Load expression data ----------
expr_raw <- read.table(
  "C:/Users/fathi/Downloads/GSE40419_LC-87_RPKM_expression.txt/GSE40419_LC-87_RPKM_expression.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ---------- Load metadata ----------
meta <- read.table(
  "C:/Users/fathi/Downloads/GSE40419_LC-87_RPKM_expression.txt/metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ---------- Define condition ----------
colnames(meta)
meta$condition <- ifelse(
  grepl("_nor", meta$sample_title, ignore.case = TRUE),
  "Normal",
  "Cancer"
)
meta$condition <- factor(meta$condition, levels = c("Normal","Cancer"))
table(meta$condition)

# ---------- Expression matrix ----------
expr <- expr_raw[, meta$sample_title]

expr <- as.matrix(expr)
mode(expr) <- "numeric"

# log2 + normalize (correct for RPKM)
expr <- log2(expr + 1)
expr <- normalizeBetweenArrays(expr)

# use gene symbols as rownames
rownames(expr) <- make.unique(expr_raw$gene)

#map technical rownames to true accessions
accession_map <- data.frame(
  row_id = rownames(expr),
  accession = expr_raw$accession,
  stringsAsFactors = FALSE
)

#LIMMA
group <- factor(meta$condition,levels = 
                  c("Normal","Tumor"))

# ---------- limma design ----------
design <- model.matrix(~ condition, data = meta)
colnames(design)

# ---------- limma fit ----------
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# ---------- Results ----------
res <- topTable(
  fit,
  coef = "conditionCancer",
  number = Inf,
  adjust.method = "BH"
)

#attach row_id
res$row_id <- rownames(res)

#merge accession info
res <- merge(res, accession_map, by = "row_id")
res$gene_symbol <- rownames(res)
colnames(res)
head(res$accession)
res$refseq <- res$accession
head(res$refseq)
res$refseq_clean <- sub("\\..*", "", res$refseq)

# ---------- biomaRt annotation ----------
mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)
annot <- getBM(
  attributes = c("refseq_mrna","hgnc_symbol"),
  filters = "refseq_mrna",
  values = unique(res$refseq_clean),
  mart = mart
)
head(annot)
#merge annotation
res <- left_join(
  res,
  annot, 
  by = c("refseq_clean" = "refseq_mrna"))
res <- res %>% filter(!is.na(hgnc_symbol))
library(dplyr)
res_unique <- res %>%
  arrange(adj.P.Val) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

sig_unique <- res_unique %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(sig_unique)


# ---------- DEG classification ----------
res$group <- "NS"
res$group[res$adj.P.Val < 0.05 & res$logFC > 1]  <- "Up"
res$group[res$adj.P.Val < 0.05 & res$logFC < -1] <- "Down"

table(res$group)

# ---------- Save DEG tables ----------
sig  <- res %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)
up   <- sig %>% filter(logFC > 1)
down <- sig %>% filter(logFC < -1)
FDR < 0.01
sig_fdr01 <- res_unique %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 1)

nrow(sig_fdr01)
sig_final <- res_unique %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2)

nrow(sig_final)


write.csv(res,  "limma_all_results.csv", row.names = FALSE)
write.csv(sig,  "DEG_significant.csv",   row.names = FALSE)
write.csv(up,   "DEG_upregulated.csv",   row.names = FALSE)
write.csv(down, "DEG_downregulated.csv", row.names = FALSE)

# ---------- PCA ----------
n_normal <- sum(meta$condition == "Normal")
n_cancer <- sum(meta$condition == "cancer")
expr_pca <- expr[rowSums(expr) > 0, ]
pca <- prcomp(t(expr_pca), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  condition = group
)
ggplot(pca_df, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle(paste0(
        "PCA – GSE40419 (Tumor n=", n_cancer,
        ", Normal n=", n_normal, ")"
    ))

# ---------- Volcano plot ----------
ggplot(res, aes(logFC, -log10(adj.P.Val), color = group)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal() +
  ggtitle("Volcano plot: Cancer vs Normal") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")

# ---------- Heatmap (Top 40 DEGs) ----------
top <- sig %>%
  arrange(adj.P.Val) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  head(40)

#extract expression using ROW_ID 
mat <- expr[top$row_id, ]

# now set readable genes names
rownames(mat) <- top$hgnc_symbol

pheatmap(
  mat,
  scale = "row",
  annotation_col = data.frame(
    condition = meta$condition,
    row.names = colnames(mat)
  ),
  show_colnames = FALSE,
  main = "Top 40 Differentially Expressed Genes"
)
res
nrow(res)
sig <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(sig)

data.frame(
  Criteria = c(
    "FDR<0.05 & |logFC|>1 (your current)",
    "FDR<0.05 & |logFC|>2",
    "FDR<0.01 & |logFC|>1",
    "FDR<0.01 & |logFC|>2"
  ),
  Gene_Count = c(
    nrow(res_unique %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)),
    nrow(res_unique %>% filter(adj.P.Val < 0.05 & abs(logFC) > 2)),
    nrow(res_unique %>% filter(adj.P.Val < 0.01 & abs(logFC) > 1)),
    nrow(res_unique %>% filter(adj.P.Val < 0.01 & abs(logFC) > 2))
  )
)

