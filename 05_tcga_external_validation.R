## =============================================================================
## External validation on TCGA-LUAD (UCSC Xena HiSeqV2: log2(RSEM norm count + 1))
## =============================================================================
## - Models are trained ONCE on all of GSE40419 and then frozen.
## - TCGA is used only for evaluation (plus label-free per-gene z-scoring).
## - Test set A: matched tumor/normal pairs. Test set B: all primary tumors vs
##   all solid-tissue normals (imbalanced, so see balanced accuracy and AP_normal).
## =============================================================================

for (p in c("limma", "glmnet", "randomForest", "pROC", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install package first: ", p)
}
suppressPackageStartupMessages({
  library(limma); library(glmnet); library(randomForest); library(pROC); library(dplyr)
})

EXPR_PATH <- "GSE40419_LC-87_RPKM_expression.txt"
PANEL     <- c("AGER", "DCN", "IL6", "LCN2", "MMP1", "MMP7", "MMP9", "SPP1")
OUT_DIR   <- "deg_outputs/external_validation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

## ---- 1. GSE40419 (training cohort), same preprocessing as 04 -----------------
raw <- read.delim(EXPR_PATH, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
ids <- make.unique(as.character(raw[[1]]))
annot_cols  <- c("accession", "chrom", "start", "end", "strand")
sample_cols <- setdiff(colnames(raw), c(colnames(raw)[1], annot_cols))
expr <- as.matrix(raw[, sample_cols]); rownames(expr) <- ids
y       <- ifelse(grepl("_nor$", sample_cols), 0, 1)
patient <- sub("_nor$", "", sample_cols)

expr <- normalizeBetweenArrays(log2(expr + 1), method = "quantile")
base_gene <- sub("\\.[0-9]+$", "", rownames(expr))
rmeans <- rowMeans(expr)
keep <- unlist(tapply(seq_along(base_gene), base_gene,
                      function(i) i[which.max(rmeans[i])]), use.names = FALSE)
expr <- expr[keep, ]; rownames(expr) <- base_gene[keep]

## ---- 2. TCGA-LUAD from UCSC Xena ---------------------------------------------
find_tcga <- function() {
  c(list.files(".", pattern = "HiSeqV2(\\.gz)?$", full.names = TRUE),
    list.files("tcga_download", pattern = "HiSeqV2(\\.gz)?$", full.names = TRUE, recursive = TRUE))
}
tcga_path <- find_tcga()
if (length(tcga_path) == 0) {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) stop("Install UCSCXenaTools first.")
  library(UCSCXenaTools)
  xe <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
    XenaFilter(filterDatasets = "TCGA.LUAD.sampleMap/HiSeqV2$")
  print(xe)
  XenaDownload(XenaQuery(xe), destdir = "tcga_download", force = FALSE)
  tcga_path <- find_tcga()
}
stopifnot("TCGA file not found; see manual download notes" = length(tcga_path) > 0)
tcga_path <- tcga_path[1]
if (grepl("\\.gz$", tcga_path) && file.size(tcga_path) < 30000000)
  stop("TCGA file looks truncated (partial download). Delete it and download again.")
cat("Using TCGA file:", tcga_path, "\n")

tc_raw <- read.delim(tcga_path, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
tc <- as.matrix(tc_raw[, -1]); rownames(tc) <- make.unique(as.character(tc_raw[[1]]))

samp_type  <- sub("^TCGA-[^-]+-[^-]+-([0-9]{2}).*$", "\\1", colnames(tc))
tc_patient <- sub("^(TCGA-[^-]+-[^-]+)-.*$", "\\1", colnames(tc))
keep_s <- samp_type %in% c("01", "11")            # 01 = primary tumor, 11 = solid normal
tc <- tc[, keep_s]; samp_type <- samp_type[keep_s]; tc_patient <- tc_patient[keep_s]
tc_y <- as.integer(samp_type == "01")
dup <- duplicated(paste(tc_patient, tc_y))         # one sample per patient per class
tc <- tc[, !dup]; tc_y <- tc_y[!dup]; tc_patient <- tc_patient[!dup]

both <- intersect(tc_patient[tc_y == 1], tc_patient[tc_y == 0])
idxA <- which(tc_patient %in% both)                # matched pairs
idxB <- seq_along(tc_y)                            # all primary tumors vs all normals
cat(sprintf("TCGA-LUAD: %d tumors, %d normals; %d matched pairs\n",
            sum(tc_y == 1), sum(tc_y == 0), length(both)))

shared <- intersect(rownames(expr), rownames(tc))
missing_panel <- setdiff(PANEL, shared)
if (length(missing_panel) > 0) stop("Panel genes missing in TCGA: ", paste(missing_panel, collapse = ", "))

## ---- 3. Helpers --------------------------------------------------------------
zmat <- function(M) sweep(sweep(M, 2, colMeans(M), "-"), 2, pmax(apply(M, 2, sd), 1e-8), "/")
inner_folds <- function(pat, k = 5) {
  u <- unique(pat); sample(rep(1:k, length.out = length(u)))[match(pat, u)]
}
fast_t <- function(E, yy) {
  a <- E[, yy == 1, drop = FALSE]; b <- E[, yy == 0, drop = FALSE]
  tt <- (rowMeans(a) - rowMeans(b)) /
    sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b))
  tt[!is.finite(tt)] <- 0
  tt
}
pr_auc <- function(yy, p) {
  o <- order(p, decreasing = TRUE); yy <- yy[o]
  sum((cumsum(yy) / seq_along(yy))[yy == 1]) / sum(yy)
}
ext_metrics <- function(yy, p, thr = 0.5) {
  pred <- as.integer(p >= thr)
  tp <- sum(pred == 1 & yy == 1); fp <- sum(pred == 1 & yy == 0)
  fn <- sum(pred == 0 & yy == 1); tn <- sum(pred == 0 & yy == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  prec <- if (tp + fp == 0) NA else tp / (tp + fp)
  f1 <- if (is.na(prec) || prec + sens == 0) NA else 2 * prec * sens / (prec + sens)
  data.frame(AUC = as.numeric(auc(roc(yy, p, levels = c(0, 1), direction = "<", quiet = TRUE))),
             AP_tumor = pr_auc(yy, p), AP_normal = pr_auc(1 - yy, 1 - p),
             Accuracy = (tp + tn) / length(yy), BalancedAcc = (sens + spec) / 2,
             Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1,
             TP = tp, FP = fp, FN = fn, TN = tn)
}
fit_models <- function(Xtr, ytr, pat) {
  d <- as.data.frame(Xtr); d$y <- ytr
  foldid <- inner_folds(pat, 5)
  list(
    Logistic     = suppressWarnings(glm(y ~ ., data = d, family = binomial)),
    LASSO        = cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1,   foldid = foldid, type.measure = "auc"),
    ElasticNet   = cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.5, foldid = foldid, type.measure = "auc"),
    RandomForest = randomForest(x = Xtr, y = factor(ytr), ntree = 500)
  )
}
predict_models <- function(f, Xte) {
  cbind(
    Logistic     = suppressWarnings(as.numeric(predict(f$Logistic, newdata = as.data.frame(Xte), type = "response"))),
    LASSO        = as.numeric(predict(f$LASSO, Xte, s = "lambda.min", type = "response")),
    ElasticNet   = as.numeric(predict(f$ElasticNet, Xte, s = "lambda.min", type = "response")),
    RandomForest = as.numeric(predict(f$RandomForest, Xte, type = "prob")[, "1"])
  )
}

## ---- 4. Train on GSE40419 (frozen), evaluate on TCGA -------------------------
tv   <- fast_t(expr[shared, ], y)                  # uses GSE40419 labels only
top8 <- names(sort(abs(tv), decreasing = TRUE))[1:8]
cat("Top-8 genes by |t| in GSE40419 (present in TCGA):", paste(top8, collapse = ", "), "\n")

gene_sets  <- list(panel8 = PANEL, top8_by_t = top8)
test_sets  <- list(matched_pairs = idxA, all_primary_vs_normal = idxB)
results <- list(); preds_out <- list()

for (gs in names(gene_sets)) {
  genes <- gene_sets[[gs]]
  Xtr  <- zmat(t(expr[genes, , drop = FALSE]))     # z-scored within GSE40419
  fits <- fit_models(Xtr, y, patient)
  for (ts in names(test_sets)) {
    idx <- test_sets[[ts]]
    Xte <- zmat(t(tc[genes, idx, drop = FALSE]))   # z-scored within the TCGA set, no labels
    pr  <- predict_models(fits, Xte)
    preds_out[[paste(gs, ts, sep = "|")]] <-
      data.frame(gene_set = gs, test_set = ts, sample = colnames(tc)[idx], y = tc_y[idx], pr)
    for (m in colnames(pr))
      results[[length(results) + 1]] <-
      cbind(gene_set = gs, test_set = ts, model = m, ext_metrics(tc_y[idx], pr[, m]))
  }
}
res_df <- do.call(rbind, results)
write.csv(res_df, file.path(OUT_DIR, "tcga_external_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, preds_out), file.path(OUT_DIR, "tcga_predictions.csv"), row.names = FALSE)

cat("\n=== TCGA-LUAD external validation ===\n")
print(res_df[, c("gene_set", "test_set", "model", "AUC", "BalancedAcc", "Sensitivity",
                 "Specificity", "AP_tumor", "AP_normal")], digits = 3, row.names = FALSE)

## ---- 5. ROC plot (panel8) ----------------------------------------------------
models <- c("LASSO", "ElasticNet", "Logistic", "RandomForest")
cols <- c("red", "blue", "black", "darkgreen")
png(file.path(OUT_DIR, "tcga_roc_panel8.png"), width = 2400, height = 1200, res = 250)
par(mfrow = c(1, 2))
for (ts in names(test_sets)) {
  d <- preds_out[[paste("panel8", ts, sep = "|")]]
  plot(roc(d$y, d[[models[1]]], levels = c(0, 1), direction = "<", quiet = TRUE),
       col = cols[1], main = paste("TCGA-LUAD,", ts))
  for (i in 2:length(models))
    plot(roc(d$y, d[[models[i]]], levels = c(0, 1), direction = "<", quiet = TRUE),
         col = cols[i], add = TRUE)
  legend("bottomright", legend = models, col = cols, lwd = 2, cex = 0.7)
}
dev.off()

## ---- 6. Do the panel genes move the same way in TCGA? (matched pairs) -------
conc <- do.call(rbind, lapply(PANEL, function(g) {
  tum <- sapply(both, function(p) tc[g, which(tc_patient == p & tc_y == 1)[1]])
  nor <- sapply(both, function(p) tc[g, which(tc_patient == p & tc_y == 0)[1]])
  data.frame(gene = g,
             log2FC_GSE40419 = mean(expr[g, y == 1]) - mean(expr[g, y == 0]),
             log2FC_TCGA_paired = mean(tum - nor),
             paired_t_p = t.test(tum, nor, paired = TRUE)$p.value)
}))
conc$same_direction <- sign(conc$log2FC_GSE40419) == sign(conc$log2FC_TCGA_paired)
conc$paired_t_FDR <- p.adjust(conc$paired_t_p, method = "BH")
write.csv(conc, file.path(OUT_DIR, "panel_gene_direction_concordance.csv"), row.names = FALSE)
cat("\n=== Panel gene direction, GSE40419 vs TCGA (matched pairs) ===\n")
print(conc, digits = 3, row.names = FALSE)

cat("\nDone. Outputs in:", normalizePath(OUT_DIR), "\n")

source("05_tcga_external_validation.R")
