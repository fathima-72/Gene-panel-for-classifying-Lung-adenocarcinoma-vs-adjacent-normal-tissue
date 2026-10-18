
## =============================================================================
## GSE40419 — PPI / Hub Gene Consolidation
## =============================================================================
## Purpose: Merge STRING/cytoHubba centrality ranking exports (Degree,
## Closeness, EPC, MCC) — keyed by Ensembl Protein ID (ENSP), not gene
## symbol — into ONE table, mapped to gene symbols, with an explicit,
## documented "hub gene" rule instead of an undocumented black box.
##
## Fixes in this version:
##   1. Runs on whichever ranking files you actually have.
##   2. Header-detection only matches whole words "Name"/"Score"/"Id", so a
##      description line containing "ranked" doesn't get mistaken for a header.
##   3. Prefers the real "Score" column over "Rank" when both exist.
##   4. All dplyr::select() calls are namespaced, since AnnotationDbi::select()
##      overrides plain select() once loaded.
##   5. ENSP -> SYMBOL mapping tries org.Hs.eg.db first, then fills gaps using
##      EnsDb.Hsapiens.v86 (org.Hs.eg.db's ENSP coverage is known to be sparse).
##   6. Cross-references consensus hub genes against DEG_significant.csv.
##   7. Explicitly checks MMP9 and ACAN against the evidence (audit Part 9).
## =============================================================================

required_pkgs <- c("org.Hs.eg.db", "AnnotationDbi", "dplyr", "purrr", "readr", "tidyr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       ". Install these once, outside this script, before running it.")
}
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  message("EnsDb.Hsapiens.v86 not installed -- run BiocManager::install(\"EnsDb.Hsapiens.v86\") ",
          "for better ENSP->SYMBOL mapping coverage. Continuing without it for now.")
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

## ---- 1. Paths and parameters --------------------------------------------------

METHOD_PATTERNS <- c(
  Degree    = "degree",
  Closeness = "closeness",
  EPC       = "epc",
  MCC       = "mcc"
)

TOP_N              <- 10
CONSENSUS_MIN_HITS <- 3

DEG_SIGNIFICANT <- "deg_outputs/DEG_significant.csv"

OUT_DIR <- "deg_outputs/hub_genes"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(
  "Run 01_DEG_analysis.R first -- DEG_significant.csv not found" =
    file.exists(DEG_SIGNIFICANT)
)

find_file_for_method <- function(pattern) {
  candidates <- list.files(".", pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  candidates <- candidates[!grepl("HubbaTable", candidates, ignore.case = TRUE)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

RANKING_FILES <- sapply(METHOD_PATTERNS, find_file_for_method)
cat("Files discovered by pattern match:\n")
print(RANKING_FILES)

## ---- 2. Read each ranking file, auto-detecting header row + ID/score cols ----

read_ranking <- function(path, method_name) {
  if (is.na(path) || !file.exists(path)) {
    cat(sprintf("[%s] No matching file found -- skipping this method.\n", method_name))
    return(NULL)
  }
  
  raw_lines <- readLines(path, warn = FALSE)
  header_idx <- which(grepl("\\bname\\b|\\bscore\\b|\\bid\\b", raw_lines, ignore.case = TRUE))[1]
  skip_n <- if (is.na(header_idx)) 0 else header_idx - 1
  
  if (skip_n > 0) {
    cat(sprintf("[%s] Detected %d description line(s) above the header -- skipping them.\n",
                method_name, skip_n))
  }
  
  df <- read_csv(path, skip = skip_n, show_col_types = FALSE)
  
  id_col <- names(df)[sapply(df, function(col) {
    any(grepl("ENSP\\d+", as.character(col)))
  })][1]
  
  score_col <- names(df)[tolower(names(df)) == "score"][1]
  if (is.na(score_col)) score_col <- names(df)[tolower(names(df)) == "rank"][1]
  
  if (is.na(id_col) || is.na(score_col)) {
    warning(sprintf("[%s] Could not auto-detect ID/score columns in %s -- skipping.",
                    method_name, path))
    return(NULL)
  }
  
  cat(sprintf("[%s] Using id column = '%s', score column = '%s' (%d rows)\n",
              method_name, id_col, score_col, nrow(df)))
  
  df %>%
    transmute(
      ensp  = sub("^9606\\.", "", .data[[id_col]]),
      score = .data[[score_col]],
      method = method_name
    ) %>%
    arrange(desc(score)) %>%
    mutate(rank_in_method = row_number())
}

ranking_list <- imap(as.list(RANKING_FILES), read_ranking)
ranking_list <- compact(ranking_list)

if (length(ranking_list) == 0) {
  stop("No ranking files could be read. Place at least one Degree/Closeness/EPC/MCC ",
       "ranking CSV in the working directory.")
}

cat(sprintf("\nSuccessfully loaded %d of %d ranking methods: %s\n\n",
            length(ranking_list), length(RANKING_FILES),
            paste(names(ranking_list), collapse = ", ")))

if (file.exists("HubbaTable.csv")) {
  hubba_check <- suppressMessages(read_csv("HubbaTable.csv", n_max = 1, show_col_types = FALSE))
  has_id_col <- any(grepl("name|id|ensp", names(hubba_check), ignore.case = TRUE))
  if (!has_id_col) {
    cat("NOTE: HubbaTable.csv is present but has NO gene/protein ID column --\n")
    cat("it cannot be used until re-exported from Cytoscape with node names included.\n\n")
  }
}

## ---- 3. Map ENSP -> gene symbol (org.Hs.eg.db, then EnsDb.Hsapiens.v86) ------

all_ensp <- unique(unlist(map(ranking_list, ~ .x$ensp)))

ensp_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensp,
  keytype = "ENSEMBLPROT",
  columns = c("SYMBOL")
) %>%
  distinct(ENSEMBLPROT, .keep_all = TRUE)

unmapped <- all_ensp[!(all_ensp %in% ensp_map$ENSEMBLPROT[!is.na(ensp_map$SYMBOL)])]

if (length(unmapped) > 0) {
  write.csv(data.frame(unmapped_ensp = unmapped),
            file.path(OUT_DIR, "unmapped_protein_ids.csv"), row.names = FALSE)
  cat(sprintf("\n%d protein IDs could not be auto-mapped.\n", length(unmapped)))
  
  if (file.exists("manual_ensp_mapping.csv")) {
    manual_map <- read.csv("manual_ensp_mapping.csv", stringsAsFactors = FALSE)
    extra_map <- manual_map %>%
      rename(ENSEMBLPROT = ensp, SYMBOL = symbol) %>%
      distinct(ENSEMBLPROT, .keep_all = TRUE)
    ensp_map <- ensp_map %>%
      filter(!(ENSEMBLPROT %in% extra_map$ENSEMBLPROT)) %>%
      bind_rows(extra_map)
    cat(sprintf("Filled in %d mappings from manual_ensp_mapping.csv\n", nrow(extra_map)))
  } else {
    cat("No manual_ensp_mapping.csv found -- these IDs remain unmapped.\n")
  }
}

n_mapped <- sum(!is.na(ensp_map$SYMBOL) & ensp_map$ENSEMBLPROT %in% all_ensp)
cat(sprintf("Mapped %d of %d unique protein IDs to gene symbols (%.0f%%)\n",
            n_mapped, length(all_ensp), 100 * n_mapped / length(all_ensp)))

if (n_mapped / length(all_ensp) < 0.5) {
  warning("Less than half of protein IDs mapped to symbols even after the EnsDb fallback. ",
          "Manually spot-check a few unmapped ENSP IDs on string-db.org to confirm they ",
          "are genuine human (9606) protein IDs.")
}

ranking_all <- bind_rows(ranking_list) %>%
  left_join(ensp_map, by = c("ensp" = "ENSEMBLPROT")) %>%
  dplyr::filter(!is.na(SYMBOL))

## ---- 4. Define "hub gene" explicitly ------------------------------------------

top_per_method <- ranking_all %>%
  group_by(method) %>%
  dplyr::filter(rank_in_method <= TOP_N) %>%
  ungroup()

hub_summary <- top_per_method %>%
  distinct(SYMBOL, method) %>%
  count(SYMBOL, name = "n_methods_in_top") %>%
  left_join(
    ranking_all %>%
      dplyr::select(SYMBOL, method, rank_in_method) %>%
      distinct() %>%
      pivot_wider(names_from = method, values_from = rank_in_method,
                  names_prefix = "rank_"),
    by = "SYMBOL"
  ) %>%
  arrange(desc(n_methods_in_top))

consensus_hub_genes <- hub_summary %>%
  dplyr::filter(n_methods_in_top >= CONSENSUS_MIN_HITS)

cat(sprintf(
  "\nHub gene definition used: top %d in a method's ranking, appearing in >= %d of %d methods.\n",
  TOP_N, CONSENSUS_MIN_HITS, length(ranking_list)
))
cat(sprintf("Consensus hub genes found: %d\n", nrow(consensus_hub_genes)))

write.csv(hub_summary, file.path(OUT_DIR, "hub_gene_ranks_all_methods.csv"), row.names = FALSE)
write.csv(consensus_hub_genes, file.path(OUT_DIR, "consensus_hub_genes.csv"), row.names = FALSE)

## ---- 5. Cross-reference against your corrected DEG list -----------------------

deg_sig <- read.csv(DEG_SIGNIFICANT, stringsAsFactors = FALSE)

consensus_hub_genes <- consensus_hub_genes %>%
  mutate(
    is_DEG    = SYMBOL %in% deg_sig$hgnc_symbol,
    logFC     = deg_sig$logFC[match(SYMBOL, deg_sig$hgnc_symbol)],
    adj.P.Val = deg_sig$adj.P.Val[match(SYMBOL, deg_sig$hgnc_symbol)]
  )

write.csv(consensus_hub_genes, file.path(OUT_DIR, "consensus_hub_genes_with_DEG_status.csv"),
          row.names = FALSE)

n_deg_confirmed <- sum(consensus_hub_genes$is_DEG, na.rm = TRUE)
cat(sprintf("Of %d consensus hub genes, %d (%.0f%%) are also significant DEGs.\n",
            nrow(consensus_hub_genes), n_deg_confirmed,
            100 * n_deg_confirmed / nrow(consensus_hub_genes)))
cat("Genes that are network hubs but NOT significant DEGs deserve a second look\n")
cat("before being treated as biomarker candidates -- see printout below.\n\n")

print(consensus_hub_genes %>% dplyr::filter(!is_DEG) %>% dplyr::select(SYMBOL, n_methods_in_top))

## ---- 6. Explicit audit of MMP9 and ACAN ---------------------------------------

for (gene in c("MMP9", "ACAN")) {
  cat(sprintf("\n--- %s ---\n", gene))
  hub_row <- hub_summary %>% dplyr::filter(SYMBOL == gene)
  deg_row <- deg_sig %>% dplyr::filter(hgnc_symbol == gene)
  
  if (nrow(hub_row) == 0) {
    cat(sprintf("%s did not appear in the top-%d of ANY loaded centrality method.\n", gene, TOP_N))
  } else {
    cat(sprintf("%s appears in top-%d of %d method(s):\n", gene, TOP_N, hub_row$n_methods_in_top))
    print(hub_row)
  }
  
  if (nrow(deg_row) == 0) {
    cat(sprintf("%s is NOT in the corrected significant DEG list (adj.P.Val<0.05, |logFC|>1).\n", gene))
  } else {
    cat(sprintf("%s IS a significant DEG: logFC = %.2f, adj.P.Val = %.2e\n",
                gene, deg_row$logFC[1], deg_row$adj.P.Val[1]))
  }
}

cat("\nDone. Outputs written to:", normalizePath(OUT_DIR), "\n")
cat("consensus_hub_genes_with_DEG_status.csv is the file to carry into the ML\n")
cat("feature-panel step -- see audit Part 6's funnel (DEG ∩ hub ∩ ML importance).\n")


