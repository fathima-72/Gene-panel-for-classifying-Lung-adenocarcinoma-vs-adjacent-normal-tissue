## =============================================================================
## GSE40419 — Leakage-aware ML classifier (Cancer vs Normal) on the 8-gene panel
## =============================================================================
## - Folds are grouped by PATIENT (LC_C1 and LC_C1_nor are the same patient).
## - Scaling is fit on the training fold only.
## - LASSO / Elastic Net tune lambda by an inner, patient-grouped CV.
## - "Logistic_infold_top8" picks its 8 genes by t-test INSIDE each training
##   fold from all genes, so it shows how optimistic the fixed-panel numbers are.
## - Feature importance at the end is descriptive (fit on all data), not a
##   performance estimate.
## =============================================================================

for (p in c("limma", "glmnet", "randomForest", "pROC", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install package first: ", p)
}
suppressPackageStartupMessages({
  library(limma); library(glmnet); library(randomForest); library(pROC); library(dplyr)
})

EXPR_PATH <- "GSE40419_LC-87_RPKM_expression.txt"
PANEL     <- c("AGER", "DCN", "IL6", "LCN2", "MMP1", "MMP7", "MMP9", "SPP1")
K         <- 5
N_REPEATS <- 3
OUT_DIR   <- "deg_outputs/ml"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and preprocess (same steps as 01_DEG_analysis.R) ----------------
raw <- read.delim(EXPR_PATH, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
ids <- make.unique(as.character(raw[[1]]))
annot_cols  <- c("accession", "chrom", "start", "end", "strand")
sample_cols <- setdiff(colnames(raw), c(colnames(raw)[1], annot_cols))
expr <- as.matrix(raw[, sample_cols]); rownames(expr) <- ids

y       <- ifelse(grepl("_nor$", sample_cols), 0, 1)   # 1 = Cancer
patient <- sub("_nor$", "", sample_cols)
cat(sprintf("%d samples (%d Cancer, %d Normal), %d patients\n",
            length(y), sum(y == 1), sum(y == 0), length(unique(patient))))

expr <- normalizeBetweenArrays(log2(expr + 1), method = "quantile")  # label-free
base_gene <- sub("\\.[0-9]+$", "", rownames(expr))
rmeans <- rowMeans(expr)
keep <- unlist(tapply(seq_along(base_gene), base_gene,
                      function(i) i[which.max(rmeans[i])]), use.names = FALSE)
expr <- expr[keep, ]; rownames(expr) <- base_gene[keep]
stopifnot(all(PANEL %in% rownames(expr)))
X <- t(expr[PANEL, ])                                  # samples x panel genes

## ---- 2. Helpers --------------------------------------------------------------
scale_fit   <- function(M) list(m = colMeans(M), s = pmax(apply(M, 2, sd), 1e-8))
apply_scale <- function(M, sc) sweep(sweep(M, 2, sc$m, "-"), 2, sc$s, "/")

inner_folds <- function(pat, k = 5) {
  u <- unique(pat)
  sample(rep(1:k, length.out = length(u)))[match(pat, u)]
}
fit_glm <- function(Xtr, ytr, Xte) {
  d <- as.data.frame(Xtr); d$y <- ytr
  fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
  suppressWarnings(as.numeric(predict(fit, newdata = as.data.frame(Xte), type = "response")))
}
fit_glmnet <- function(Xtr, ytr, Xte, alpha, foldid) {
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = alpha,
                  foldid = foldid, type.measure = "auc")
  as.numeric(predict(cv, Xte, s = "lambda.min", type = "response"))
}
fit_rf <- function(Xtr, ytr, Xte) {
  rf <- randomForest(x = Xtr, y = factor(ytr), ntree = 500)
  predict(rf, Xte, type = "prob")[, "1"]
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
  sum((cumsum(yy) / seq_along(yy))[yy == 1]) / sum(yy)   # average precision
}
metrics <- function(yy, p, thr = 0.5) {
  pred <- as.integer(p >= thr)
  tp <- sum(pred == 1 & yy == 1); fp <- sum(pred == 1 & yy == 0)
  fn <- sum(pred == 0 & yy == 1); tn <- sum(pred == 0 & yy == 0)
  prec <- if (tp + fp == 0) NA else tp / (tp + fp)
  rec  <- tp / (tp + fn)
  f1   <- if (is.na(prec) || prec + rec == 0) NA else 2 * prec * rec / (prec + rec)
  data.frame(AUC = as.numeric(auc(roc(yy, p, levels = c(0, 1), direction = "<", quiet = TRUE))),
             PR_AUC = pr_auc(yy, p), Accuracy = (tp + tn) / length(yy),
             Precision = prec, Recall = rec, F1 = f1, TP = tp, FP = fp, FN = fn, TN = tn)
}

## ---- 3. Patient-grouped cross-validation -------------------------------------
model_names <- c("Logistic", "LASSO", "ElasticNet", "RandomForest", "Logistic_infold_top8")

run_repeat <- function(rep_i) {
  set.seed(1000 + rep_i)
  u <- unique(patient)
  fold <- sample(rep(1:K, length.out = length(u)))[match(patient, u)]
  oof <- matrix(NA_real_, nrow = length(y), ncol = length(model_names),
                dimnames = list(sample_cols, model_names))
  chosen <- list()
  for (k in 1:K) {
    tr <- which(fold != k); te <- which(fold == k)
    sc  <- scale_fit(X[tr, , drop = FALSE])
    Xtr <- apply_scale(X[tr, , drop = FALSE], sc)
    Xte <- apply_scale(X[te, , drop = FALSE], sc)
    inner <- inner_folds(patient[tr])
    oof[te, "Logistic"]     <- fit_glm(Xtr, y[tr], Xte)
    oof[te, "LASSO"]        <- fit_glmnet(Xtr, y[tr], Xte, 1, inner)
    oof[te, "ElasticNet"]   <- fit_glmnet(Xtr, y[tr], Xte, 0.5, inner)
    oof[te, "RandomForest"] <- fit_rf(Xtr, y[tr], Xte)
    
    # leakage-free comparison: choose the top 8 genes inside the training fold
    tv  <- fast_t(expr[, tr], y[tr])
    top <- names(sort(abs(tv), decreasing = TRUE))[seq_along(PANEL)]
    chosen[[k]] <- top
    X2  <- t(expr[top, , drop = FALSE])
    sc2 <- scale_fit(X2[tr, , drop = FALSE])
    oof[te, "Logistic_infold_top8"] <- fit_glm(apply_scale(X2[tr, , drop = FALSE], sc2), y[tr],
                                               apply_scale(X2[te, , drop = FALSE], sc2))
  }
  list(fold = fold, oof = oof, chosen = chosen)
}

res <- lapply(seq_len(N_REPEATS), run_repeat)

## ---- 4. Metrics --------------------------------------------------------------
all_metrics <- do.call(rbind, lapply(seq_along(res), function(i)
  do.call(rbind, lapply(model_names, function(m)
    cbind(repeat_id = i, model = m, metrics(y, res[[i]]$oof[, m]))))))

summ <- all_metrics %>%
  group_by(model) %>%
  summarise(across(c(AUC, PR_AUC, Accuracy, Precision, Recall, F1),
                   list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))),
            .groups = "drop")

write.csv(all_metrics, file.path(OUT_DIR, "cv_metrics_per_repeat.csv"), row.names = FALSE)
write.csv(summ, file.path(OUT_DIR, "cv_metrics_summary.csv"), row.names = FALSE)

oof_df <- data.frame(sample = sample_cols, patient = patient, y = y,
                     fold_repeat1 = res[[1]]$fold, res[[1]]$oof)
write.csv(oof_df, file.path(OUT_DIR, "oof_predictions_repeat1.csv"), row.names = FALSE)

cat("\n=== Patient-grouped CV (mean over repeats) ===\n")
print(as.data.frame(summ), digits = 3)
cat("\nConfusion matrix, LASSO, repeat 1 (rows = truth, cols = predicted at 0.5):\n")
print(table(truth = y, predicted = as.integer(res[[1]]$oof[, "LASSO"] >= 0.5)))

cat("\nGenes chosen by the in-fold t-test (count out of", N_REPEATS * K, "folds):\n")
sel <- sort(table(unlist(lapply(res, function(r) unlist(r$chosen)))), decreasing = TRUE)
print(head(sel, 20))
cat("Overlap with your fixed panel:", paste(intersect(names(sel), PANEL), collapse = ", "), "\n")

## ---- 5. ROC plot (repeat 1) --------------------------------------------------
cols <- c("black", "red", "blue", "darkgreen", "orange")
png(file.path(OUT_DIR, "roc_curves_repeat1.png"), width = 1600, height = 1400, res = 250)
plot(roc(y, res[[1]]$oof[, 1], levels = c(0, 1), direction = "<", quiet = TRUE),
     col = cols[1], main = "Out-of-fold ROC (patient-grouped 5-fold CV)")
for (i in 2:length(model_names))
  plot(roc(y, res[[1]]$oof[, i], levels = c(0, 1), direction = "<", quiet = TRUE),
       col = cols[i], add = TRUE)
legend("bottomright", legend = model_names, col = cols, lwd = 2, cex = 0.8)
dev.off()

## ---- 6. Descriptive feature importance (all data; NOT a performance estimate) --
Xs <- apply_scale(X, scale_fit(X))
d  <- as.data.frame(Xs); d$y <- y
glm_fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
cvfit   <- cv.glmnet(Xs, y, family = "binomial", alpha = 1,
                     foldid = inner_folds(patient), type.measure = "auc")
rf_fit  <- randomForest(x = Xs, y = factor(y), ntree = 500, importance = TRUE)
imp <- data.frame(
  gene = PANEL,
  logistic_abs_coef = abs(coef(glm_fit)[PANEL]),
  lasso_abs_coef    = abs(as.numeric(coef(cvfit, s = "lambda.min"))[-1]),
  rf_gini           = as.numeric(importance(rf_fit)[PANEL, "MeanDecreaseGini"])
)
imp$mean_rank <- rowMeans(cbind(rank(-imp$logistic_abs_coef), rank(-imp$lasso_abs_coef), rank(-imp$rf_gini)))
imp <- imp[order(imp$mean_rank), ]
write.csv(imp, file.path(OUT_DIR, "feature_importance_descriptive.csv"), row.names = FALSE)
cat("\n=== Feature importance (descriptive) ===\n"); print(imp, digits = 3, row.names = FALSE)

cat("\nDone. Outputs in:", normalizePath(OUT_DIR), "\n")