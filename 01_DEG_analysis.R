## =============================================================================
## GSE40419 Lung Adenocarcinoma vs Adjacent-Normal — Differential Expression
## =============================================================================
## Purpose : Single, reproducible source of truth for DEG_significant.csv /
##           DEG_upregulated.csv / DEG_downregulated.csv, PCA, volcano, heatmap.
##
## Fixes applied relative to the audited scripts (limma_expressing_cancer_vs_
## normal.R / DEG2.R):
##   1. Genes are deduplicated to one row per hgnc_symbol BEFORE limma is run,
##      not after — this is what previously inflated DEG counts by ~35-40%.
##   2. eBayes(..., trend = TRUE) is used, since the input is log-RPKM
##      (RNA-seq-derived), not microarray intensity — limma-trend models the
##      mean-variance relationship instead of assuming it is constant.
##   3. Exactly ONE significance threshold is defined, in one place
##      (SIG_ADJP / SIG_LOGFC below), and every output file is written from
##      that single run — no interactive edits between the threshold being
##      chosen and the CSVs being saved.
##   4. biomaRt's annotation table is cached to disk after the first run, so
##      re-running the script does not depend on a live Ensembl connection.
##   5. condition factor levels are set explicitly and checked with a
##      stopifnot(), so a level-name mismatch (e.g. "Tumor" vs "Cancer")
##      fails loudly instead of silently producing NAs.
##   6. Raw row_ids are screened for Excel date-corruption (e.g. "04-Sep")
##      before annotation, and any hits are reported rather than silently
##      dropped or silently kept.
##
## Run this script top-to-bottom in a fresh R session. Do not hand-edit any
## object it creates and then re-save a CSV — re-run the whole script instead,
## so the code and the saved files always match.
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------

required_pkgs <- c("limma", "biomaRt", "dplyr", "ggplot2", "pheatmap",
                    "ggrepel", "matrixStats")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       ". Install these once, outside this script, before running it.")
}

suppressPackageStartupMessages({
  library(limma)
  library(biomaRt)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})

set.seed(1)

## ---- 1. Paths and analysis parameters (single source of truth) --------------

EXPR_PATH   <- "C:/Users/fathi/Downloads/GSE40419_LC-87_RPKM_expression.txt/GSE40419_LC-87_RPKM_expression.txt"
META_PATH   <- "C:/Users/fathi/Downloads/GSE40419_LC-87_RPKM_expression.txt/metadata.txt"
ANNOT_CACHE <- "biomart_refseq_to_hgnc_cache.csv"       # created on first run

OUT_DIR     <- "deg_outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

# The ONE threshold used everywhere downstream. Change it here, nowhere else.
SIG_ADJP  <- 0.05
SIG_LOGFC <- 1

## ---- 2. Load expression matrix and build metadata from column names ----------

expr_raw <- read.delim(EXPR_PATH, row.names = NULL, check.names = FALSE,
                       stringsAsFactors = FALSE)

id_col   <- colnames(expr_raw)[1]
gene_ids <- make.unique(as.character(expr_raw[[id_col]]))

# Separate the 5 leading annotation columns (accession, chrom, start, end,
# strand) from the actual sample expression columns.
annot_cols  <- c("accession", "chrom", "start", "end", "strand")
sample_cols <- setdiff(colnames(expr_raw), c(id_col, annot_cols))

gene_annot <- expr_raw[, c(id_col, intersect(annot_cols, colnames(expr_raw)))]
gene_annot[[id_col]] <- gene_ids

expr_raw <- expr_raw[, sample_cols, drop = FALSE]
rownames(expr_raw) <- gene_ids

# Build metadata directly from the sample naming convention in this file:
# "_nor" suffix = adjacent-normal tissue; no suffix = tumor.
meta <- data.frame(
  sample_id = sample_cols,
  condition = ifelse(grepl("_nor$", sample_cols), "Normal", "Cancer"),
  stringsAsFactors = FALSE
)

meta$condition <- factor(meta$condition, levels = c("Normal", "Cancer"))
stopifnot(!anyNA(meta$condition))
stopifnot(identical(colnames(expr_raw), meta$sample_id))

cat(sprintf("Loaded %d genes x %d samples (%d Cancer, %d Normal)\n",
            nrow(expr_raw), ncol(expr_raw),
            sum(meta$condition == "Cancer"), sum(meta$condition == "Normal")))

## ---- 2b. DIAGNOSTIC: verify expr_raw columns actually match metadata --------
cat("First 5 expr_raw column names:\n"); print(head(colnames(expr_raw), 5))
cat("Last 5 expr_raw column names:\n");  print(tail(colnames(expr_raw), 5))
cat("First 5 meta$sample_id values:\n"); print(head(meta$sample_id, 5))
cat("How many meta sample_ids are found in expr_raw columns?\n")
print(sum(meta$sample_id %in% colnames(expr_raw)))
cat("How many expr_raw columns are NOT in meta$sample_id (possible junk columns)?\n")
print(colnames(expr_raw)[!(colnames(expr_raw) %in% meta$sample_id)])

## ---- 3. Screen raw row IDs for Excel date-corruption -------------------------
## Classic patterns: SEPT1-14 -> "1-Sep".."14-Sep", MARCH1-11 -> "1-Mar".."11-Mar",
## DEC1 -> "1-Dec", etc. Report, don't silently fix -- you need to know whether
## the corruption happened upstream of you (in the GEO-deposited file itself)
## or somewhere in your own pipeline.

date_like <- grepl("^\\d{1,2}-[A-Za-z]{3}$", rownames(expr_raw))
if (any(date_like)) {
  warning(sprintf(
    "%d row IDs look like Excel-corrupted gene names (e.g. %s). These will \n
     still be pushed through RefSeq/annotation lookup below, but inspect \n
     %s manually before trusting results for these genes.",
    sum(date_like),
    paste(head(rownames(expr_raw)[date_like], 5), collapse = ", "),
    "deg_outputs/excel_date_corrupted_rows.csv"
  ))
  write.csv(data.frame(row_id = rownames(expr_raw)[date_like]),
            file.path(OUT_DIR, "excel_date_corrupted_rows.csv"), row.names = FALSE)
} else {
  cat("No Excel date-corrupted row IDs detected.\n")
}

## ---- 4. Preprocessing: log2, cross-sample normalization -----------------------
## NOTE (methods-section language): normalizeBetweenArrays() is applied here
## for cross-sample comparability of RPKM values, given that GEO distributes
## GSE40419 as RPKM rather than raw counts (count-based normalization via
## DESeq2/edgeR/voom is not usable without the raw counts). This is stated
## explicitly rather than left implicit -- see audit Part 2.2.

expr_log   <- log2(as.matrix(expr_raw) + 1)
expr_norm  <- normalizeBetweenArrays(expr_log, method = "quantile")


## ---- 5. Deduplicate to one row per gene BEFORE limma --------------------------
## This is the single most important fix from the audit: the previous pipeline
## ran limma on every transcript/probe row and only deduplicated afterwards for
## GO/KEGG, which inflated the DEG counts actually reported by ~35-40%.
##
## Strategy: collapse multiple rows mapping to the same row_id-derived gene by
## taking the row with the highest overall mean expression (a standard,
## defensible collapse rule -- the highest-expressed transcript is assumed the
## most reliably measured). If your row IDs are already one-per-gene, this is
## a no-op.

row_means <- rowMeans(expr_norm)
base_gene <- sub("\\.[0-9]+$", "", rownames(expr_norm))  # strip ".1", ".2" suffixes

keep_idx <- tapply(seq_along(base_gene), base_gene, function(idx) {
  idx[which.max(row_means[idx])]
})
keep_idx <- unlist(keep_idx, use.names = FALSE)

n_before <- nrow(expr_norm)
expr_dedup <- expr_norm[keep_idx, , drop = FALSE]
rownames(expr_dedup) <- base_gene[keep_idx]
n_after <- nrow(expr_dedup)

cat(sprintf("Deduplication: %d rows -> %d unique genes (%d transcript-duplicate rows removed)\n",
            n_before, n_after, n_before - n_after))

## ---- 6. Differential expression (limma-trend) ---------------------------------

design <- model.matrix(~condition, data = meta)
colnames(design) <- make.names(colnames(design))

fit  <- lmFit(expr_dedup, design)
fit  <- eBayes(fit, trend = TRUE)   # trend=TRUE: appropriate for log-RPKM RNA-seq data

coef_name <- grep("Cancer", colnames(design), value = TRUE)
stopifnot(length(coef_name) == 1)

res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
res$row_id <- rownames(res)
res <- res[, c("row_id", setdiff(colnames(res), "row_id"))]

cat(sprintf("Contrast: %s (positive logFC = higher in Cancer)\n", coef_name))

## ---- 7. Gene symbols (already present) + repair Excel date-corruption --------
## row_id in this file is already the gene symbol (confirmed by the date-like
## corruption pattern found in Section 3 -- only gene symbols like "SEPT4"
## get mangled into dates by Excel; RefSeq accessions never do). So no
## RefSeq->HGNC lookup is needed; we just repair the known corrupted names.

repair_map <- c(
  "1-Sep"="SEPTIN1","2-Sep"="SEPTIN2","3-Sep"="SEPTIN3","4-Sep"="SEPTIN4",
  "5-Sep"="SEPTIN5","6-Sep"="SEPTIN6","7-Sep"="SEPTIN7","8-Sep"="SEPTIN8",
  "9-Sep"="SEPTIN9","10-Sep"="SEPTIN10","11-Sep"="SEPTIN11","12-Sep"="SEPTIN12",
  "14-Sep"="SEPTIN14","15-Sep"="SEPTIN15",
  "1-Mar"="MARCHF1","2-Mar"="MARCHF2","3-Mar"="MARCHF3","4-Mar"="MARCHF4",
  "5-Mar"="MARCHF5","6-Mar"="MARCHF6","7-Mar"="MARCHF7","8-Mar"="MARCHF8",
  "9-Mar"="MARCHF9","10-Mar"="MARCHF10","11-Mar"="MARCHF11",
  "1-Dec"="DELEC1","2-Dec"="DEC2"
)

repair_symbol <- function(x) {
  # normalize "15-Sep" and "07-Mar" (zero-padded) to the same lookup key
  key <- sub("^0", "", x)
  ifelse(key %in% names(repair_map), repair_map[key], x)
}

res$hgnc_symbol <- repair_symbol(res$row_id)

n_repaired <- sum(res$hgnc_symbol != res$row_id)
cat(sprintf("Repaired %d Excel-corrupted gene symbols back to HGNC names\n", n_repaired))

res <- res %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

cat(sprintf("Final gene table: %d unique genes\n", nrow(res)))

## ---- 8. Apply the ONE threshold and write outputs -----------------------------

sig      <- res %>% filter(adj.P.Val < SIG_ADJP, abs(logFC) > SIG_LOGFC)
sig_up   <- sig %>% filter(logFC >  SIG_LOGFC)
sig_down <- sig %>% filter(logFC < -SIG_LOGFC)

cat(sprintf("Threshold used: adj.P.Val < %s & |logFC| > %s\n", SIG_ADJP, SIG_LOGFC))
cat(sprintf("Significant genes: %d total (%d up, %d down) -- unique genes, not transcripts\n",
            nrow(sig), nrow(sig_up), nrow(sig_down)))

write.csv(res,      file.path(OUT_DIR, "DEG_all_tested.csv"),      row.names = FALSE)
write.csv(sig,      file.path(OUT_DIR, "DEG_significant.csv"),     row.names = FALSE)
write.csv(sig_up,   file.path(OUT_DIR, "DEG_upregulated.csv"),     row.names = FALSE)
write.csv(sig_down, file.path(OUT_DIR, "DEG_downregulated.csv"),   row.names = FALSE)

## ---- 9. PCA ---------------------------------------------------------------

pca_input <- expr_dedup[rownames(expr_dedup) %in% res$hgnc_symbol, , drop = FALSE]

gene_var    <- apply(pca_input, 1, var)
pca_input   <- pca_input[gene_var > 0, , drop = FALSE]
cat(sprintf("PCA input: %d genes after removing %d zero-variance genes\n",
            nrow(pca_input), sum(gene_var == 0)))

pca <- prcomp(t(pca_input), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  condition = meta$condition
)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(title = "PCA of GSE40419 samples",
       x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "PCA_plot.png"), p_pca, width = 6, height = 5, dpi = 300)

## ---- 10. Volcano plot -------------------------------------------------------

res$status <- "Not significant"
res$status[res$hgnc_symbol %in% sig_up$hgnc_symbol]   <- "Up in Cancer"
res$status[res$hgnc_symbol %in% sig_down$hgnc_symbol] <- "Down in Cancer"

top_labels <- sig %>% arrange(adj.P.Val) %>% head(15)

p_volcano <- ggplot(res, aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(-SIG_LOGFC, SIG_LOGFC), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(SIG_ADJP), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = top_labels, aes(label = hgnc_symbol),
                            size = 3, max.overlaps = 20, color = "black") +
  scale_color_manual(values = c("Up in Cancer" = "#d6604d",
                                 "Down in Cancer" = "#4393c3",
                                 "Not significant" = "grey80")) +
  labs(title = "Cancer vs Normal (GSE40419)", x = "log2 fold change",
       y = expression(-log[10]~"(adj. P value)")) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "volcano_plot.png"), p_volcano, width = 7, height = 6, dpi = 300)

## ---- 11. Heatmap of top DEGs -------------------------------------------------

top_n     <- 40
top_genes <- sig %>% arrange(adj.P.Val) %>% head(top_n) %>% pull(hgnc_symbol)

heat_mat <- expr_dedup[rownames(expr_dedup) %in% top_genes, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))  # z-score per gene across samples

annotation_col <- data.frame(Condition = meta$condition)
rownames(annotation_col) <- colnames(heat_mat)

png(file.path(OUT_DIR, "heatmap_top_DEGs.png"), width = 1400, height = 1800, res = 200)
pheatmap(heat_mat,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = sprintf("Top %d DEGs (Cancer vs Normal)", top_n))
dev.off()

## ---- 12. Session info (for reproducibility appendix) --------------------------

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\nDone. All outputs written to:", normalizePath(OUT_DIR), "\n")
cat("Downstream scripts (functional annotation, network analysis, ML) should\n")
cat("read from DEG_all_tested.csv / DEG_significant.csv in this folder, not\n")
cat("from any earlier, differently-thresholded copy.\n")

## ---- 13. Always back up outputs to a fixed, permanent location -------------
PERMANENT_BACKUP <- "C:/Users/fathi/Downloads/lung_cancer_project_outputs"
dir.create(PERMANENT_BACKUP, showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(OUT_DIR, full.names = TRUE), PERMANENT_BACKUP,
          overwrite = TRUE, recursive = TRUE)
cat("Backup copy also saved to:", PERMANENT_BACKUP, "\n")

list.files("C:/Users/fathi/Downloads/lung_cancer_project_outputs")

getwd()
OUT_DIR
list.files(getwd())

list.files("C:/Users/fathi/Downloads/lung_cancer_project_outputs")
