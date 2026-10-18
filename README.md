**Network-informed gene panel for lung adenocarcinoma vs adjacent-normal tissue**
RNA-seq analysis of GSE40419 (87 tumor, 77 adjacent-normal samples): differential expression, pathway enrichment, protein-interaction network hubs, a patient-grouped machine-learning classifier, and external-validation on TCGA-LUAD.

**Key results**

-**22,428 genes tested, 1,333 significant** (485 up, 848 down; adj. P < 0.05, |log2FC| > 1). An earlier version overcounted by treating transcript rows as genes; this is fixed.
-**Pathways**: extracellular-matrix remodeling, cell adhesion, proliferation, angiogenesis, immune signaling.
-**8-gene panel** (network hubs that are also DEGs): AGER, DCN, IL6, LCN2, MMP1, MMP7, MMP9, SPP1.
-**Internal cross-validation** (folds grouped by patient): AUC 0.94 to 0.99 across four models.
-**External validation on TCGA-LUAD** (58 matched tumor/normal pairs, models frozen after training on GSE40419): AUC 0.995 to 0.999. All 8 genes change in the same direction (FDR < 0.001).
-**Honest finding**: a purely statistical top-8 (chosen by t-test inside each training fold) performs as well as the network panel. The network panel's value is interpretability, not accuracy.

## Repository layout
 
```
scripts/                      R scripts, run in order
  01_DEG_analysis.R
  02_functional_annotation.R
  03_hub_gene_consolidation.R
  04_ml_classifier.R
  05_tcga_external_validation.R
Degree.csv, Closeness.csv,   cytoHubba exports and STRING ID mapping
  MCC.csv, EPC.csv,            (inputs for script 03, kept in the project root)
  manual_ensp_mapping.csv
deg_outputs/                  results: DEG tables, enrichment, hub genes, ML, validation
docs/                         report and figures
```
## Data (not included)
 
- GSE40419: download `GSE40419_LC-87_RPKM_expression.txt` from GEO and place it in the project root.
- TCGA-LUAD: `05_tcga_external_validation.R` downloads it from UCSC Xena (about 30 MB).

## Run
 
```r
setwd("path/to/this/folder")
source("scripts/01_DEG_analysis.R")
source("scripts/02_functional_annotation.R")
source("scripts/03_hub_gene_consolidation.R")
source("scripts/04_ml_classifier.R")
source("scripts/05_tcga_external_validation.R")
```

R 4.5.2 and Bioconductor 3.22. Package versions: `deg_outputs/sessionInfo.txt`.

## Limitations (short version)
 
The panel was chosen using all discovery samples, so internal scores are optimistic and TCGA is the real test. Tumor vs adjacent normal is an easy contrast and says nothing about clinical nodule classification. Per-cohort z-scoring means a single new patient cannot be scored in isolation. See the report for the full list.
 
## License
 
MIT, see `LICENSE`.
