# Cell Type Annotation Benchmarking Framework

> This is the code-only benchmarking pipeline extracted from the thesis repository (see the
> top-level [README](../README.md) for the full project and phase-mapping context). Input
> datasets, the CellMarker2.0 reference database, and Python virtual environments are **not**
> included here — see `docs/install_notes.R` / `docs/install_notes_DL.txt` for how to
> reconstitute the environment, and the top-level README for dataset sourcing.

## Overview

This framework provides a standardized benchmarking system for evaluating 51+ cell type annotation tools using cross-validation. It supports marker-based, reference-based, classic ML, deep learning, and LLM-based approaches with a unified interface, comprehensive metrics (accuracy, F1, Cohen's Kappa), and memory/runtime tracking.

## General Workflow

### 1. Data Preparation
- Load Seurat object with single-cell RNA-seq data
- Ensure `Ground_Truth_Celltype` column exists in metadata
- Data should be preprocessed (normalized, scaled as needed by individual tools)

### 2. Cross-Validation Setup
```r
# Create k-fold cross-validation splits
folds <- create_cv_folds(seurat_obj, k = 5, group_by = "donor_id", stratify_by = "Ground_Truth_Celltype")

# Options:
# - Grouped CV: group_by = "donor_id" (keeps donor samples together)
# - Stratified CV: group_by = NA (balances cell types across folds)
```

### 3. Tool Execution Pipeline
```r
# Each fold processes:
# 1. Split data into train/test based on fold
# 2. Find marker genes on training data: FindAllMarkers(seurat_train)
# 3. Run tool function: tool_function(seurat_train, seurat_test, markers)
# 4. Collect results and calculate metrics
```

### 4. Results Aggregation
- **Pooled analysis (Method B)**: Combine all cell predictions across folds for overall metrics
- **Fold variation (Method A)**: Calculate per-fold metrics then aggregate (mean, std, 95% CI)
- **Metrics**: Accuracy, precision, recall, F1, macro F1, Cohen's Kappa
- **Resource tracking**: Runtime (wall-clock) and peak system memory per fold

## Tool Function Structure

### Standard Function Signature
```r
run_[TOOL]_function <- function(seurat_train, seurat_test, markers) {
  # Implementation here
}
```

### Required Inputs
- **`seurat_train`**: Training Seurat object (reference data)
- **`seurat_test`**: Test Seurat object (query data to predict)
- **`markers`**: Data frame from `FindAllMarkers()` output with columns:
  - `gene`: Gene symbol
  - `cluster`: Cell type/cluster identifier
  - `avg_log2FC`: Average log2 fold change
  - `p_val_adj`: Adjusted p-value
  - `pct.1`: Percentage of cells expressing gene in cluster
  - `pct.2`: Percentage of cells expressing gene in other clusters

### Required Outputs
Must return a list with exactly these elements:
```r
return(list(
  predictions = as.character(cell_predictions),      # Vector of predicted cell types
  true_labels = seurat_test$Ground_Truth_Celltype,   # Vector of true cell types
  confidence_scores = confidence_scores,             # Numeric vector of confidence scores
  cell_ids = colnames(seurat_test)                   # Vector of cell identifiers
))
```

Tools may also include `runtime_secs` and `peak_system_memory_mb` in the return list; the framework falls back to wall-clock timing if `runtime_secs` is absent or NULL.

### Implementation Guidelines

#### 1. Library Loading
```r
if (!requireNamespace("PACKAGE", quietly = TRUE)) {
  stop("PACKAGE not available. Please install it first.")
}
library(PACKAGE)
```

#### 2. Error Handling
```r
default_return <- function() {
  return(list(
    predictions = rep("Unknown", ncol(seurat_test)),
    true_labels = seurat_test$Ground_Truth_Celltype,
    confidence_scores = rep(0, ncol(seurat_test)),
    cell_ids = colnames(seurat_test)
  ))
}
```

#### 3. Memory Tracking
All tools wrap their core algorithm in `peakRAM::peakRAM()` for memory measurement. Python-based tools (SCSA, CellTypist, scAnnotate) parse `PEAK_SYSTEM_MEMORY_MB:` and `RUNTIME_SECS:` from stdout.

#### 4. Backward Compatibility
```r
run_[TOOL] <- run_[TOOL]_function
```

## Tool Categories

### Marker-Based Tools (`marker-based/`) — 10 tools
Use marker genes to classify cells. All tools call `prepare_markers()` from `benchmarking_helpers.R` for uniform filtering — see **Marker Filtering** section below for actual thresholds. SCINA is the only exception, using inline filtering instead.

| Tool | Algorithm | Notable |
|------|-----------|---------|
| **scCATCH** | Database-driven marker matching | Cluster-level predictions; uses `Ground_Truth_Celltype` as cluster vector — potential leakage |
| **SCSA** | Python-based annotation | Requires conda env; hardcoded miniforge path; `markers` param ignored (DB-lookup tool — uses test-side clustering instead) |
| **scType** | Excel-based marker database | Sources GitHub URLs per fold (network-dependent); 2-level fallback (wrapper → direct `sctype_score`) |
| **SCINA** | Semi-supervised probabilistic | Bypasses `prepare_markers()`; inline filtering (`0.5/0.05/0.15/top 20`); `allow_unknown=FALSE`; outputs `"unknown"` not `"Unknown"` |
| **scSorter** | Marker-weight sorting | Uses `prepare_markers()`; all cells assigned (no Unknown output) |
| **CellAssign** | TensorFlow probabilistic | Uses `prepare_markers()`; TF memory not captured by peakRAM; TF compatibility issues |
| **Garnett** | Monocle3 tree-based classifier | `cluster_extend=TRUE`; fallback uses ground-truth majority vote — **data leakage** (see Known Issues) |
| **ScInfeR** | kNN-weighted inference | Injects `log_warning()` into `.GlobalEnv` with `on.exit` cleanup; requires UMAP on test data |
| **clustifyr_hyper** | Hypergeometric enrichment | Cluster-level; wide-format marker matrix; binary confidence scores only |
| **clustifyr_jaccard** | Jaccard similarity | Same as hyper with different metric; binary confidence scores only |

### Reference-Based Tools (`reference-based/`) — 8 tools
Use training data as reference. Most ignore the `markers` parameter.

| Tool | Algorithm | Uses Markers? |
|------|-----------|---------------|
| **SingleR** | DE-based reference classification | Yes (as gene filter) |
| **scmap_cell** | Cell-level kNN mapping | No |
| **scmap_cluster** | Cluster-level projection | No |
| **scibetR** | TPM-based SciBet classifier | No |
| **Seurat Transfer PCA** | Transfer anchors + PCA projection | No |
| **Seurat Transfer CCA** | Transfer anchors + CCA | No |
| **Seurat Transfer RPCA** | Transfer anchors + reciprocal PCA | No |
| **CIPR** | Logfc dot product comparison | No (recomputes internally) |

### Classic ML-Based Tools (`classic-ML-based/`) — 30+ tools
Train classical ML classifiers on expression data.

- **scPred** (19 variants): svmRadial, glm, knn, rf, xgbTree, nb, nnet, multinom, pls, earth, lda, lda2, C5.0, bagEarth, treebag, AdaBoost, LogitBoost, pda, hdda
- **scAnnotatR**: SVM-based annotation
- **CALLR**: Correlation-based
- **CHETAH**: Classification tree
- **scClassify**: Ensemble learning
- **singleCellNet**: Random forest
- **CellTypist**: Python logistic regression (conda env)
- **scAnnotate**: Python-based (conda env)
- **scAnno**: kNN-based
- **CaSTLe**: XGBoost-based
- **scID**: Marker-weighted scoring

### Deep Learning-Based Tools (`DL-based/`) — 4 tools
- **CAMLU**: Deep learning annotation
- **NeuCA**: Neural network classifier
- **scLearn**: Deep learning with feature selection
- **scPred_mxnetAdam**: scPred with mxnet Adam optimizer

### LLM-Based Tools (`LLM-based/`) — 10 tools
Use LLM reasoning over marker gene lists (marker-based, no reference needed).

- **GPTCelltype**: GPT-based cell type inference
- **CASSIA**: LLM-powered annotation
- **mLLMCelltype** (8 variants): claude_sonnet4.5, deepseek_v3.1, gemini_pro2.5, gpt5, grok4_fast, llama4_maverick, qwen3_max, base

## Helper Functions (`benchmarking_helpers.R`)

### Core Functions
| Function | Purpose |
|----------|---------|
| `create_cv_folds()` | Create k-fold CV splits (grouped or stratified via caret) |
| `get_fold_data_cached()` | Split fold data + FindAllMarkers with global caching (`fold_cache`) |
| `run_single_tool_cv()` | Execute one tool across all folds with tryCatch error handling |
| `calculate_metrics()` | Accuracy, macro F1, Cohen's Kappa, per-class precision/recall/F1 |
| `aggregate_tool_results()` | Pooled + fold-variation stats, runtime, memory for one tool |
| `aggregate_all_tools()` | Comparison table across all tools, sorted by pooled accuracy |

### Metrics Details
- **calculate_metrics**: Excludes `"Unknown"` (capital U) from per-class metrics; Cohen's Kappa adjusts for chance agreement; also computes MCC and rare-type F1
- **aggregate_tool_results**: 95% CIs calculated as ±1.96 × SE; sorted by pooled kappa (not accuracy)
- **aggregate_all_tools**: Returns comparison_table, best_tool, summary_stats (num_tools, best_accuracy, best_f1, best_kappa, best_mcc, fastest_tool)

### Usage Example
```r
source("benchmarking_helpers.R")
source("marker-based/run_scCATCH.R")

folds <- create_cv_folds(seurat_obj, k = 5)
results <- run_single_tool_cv(seurat_obj, run_scCATCH_function, folds, "scCATCH")

print(results$pooled_metrics$overall_accuracy)
print(results$pooled_metrics$cohens_kappa)
```

## File Organization

```
02-benchmarking-pipeline/
├── README.md                              # This documentation
├── run_benchmark.R                        # Main orchestration (sources all tools, runs CV)
├── benchmarking_helpers.R                 # Core framework functions
├── docs/
│   ├── marker-filtering-standardization.md  # Marker filtering standards & exceptions
│   ├── memory-and-markers-audit.md          # Memory tracking & marker usage audit
│   ├── install_notes.R                      # R package install commands
│   └── install_notes_DL.txt                 # Deep-learning tool conda environment setup
├── database-based/                        # 9 database-marker tools (Phase 4 configuration;
│   │                                       # requires the CellMarker2.0 reference DB, not included)
│   └── ...
├── marker-based/                          # 10 marker-based tools (matched-reference oracle config)
│   ├── run_scCATCH.R
│   ├── run_SCSA.R
│   ├── run_scType.R
│   ├── run_SCINA.R
│   ├── run_scSorter.R
│   ├── run_CellAssign.R
│   ├── run_Garnett.R
│   ├── run_ScInfeR.R
│   ├── run_clustifyr_hyper.R
│   └── run_clustifyr_jaccard.R
├── reference-based/                       # 8 reference-based tools
│   ├── run_SingleR.R
│   ├── run_scmap_cell.R
│   ├── run_scmap_cluster.R
│   ├── run_scibetR.R
│   ├── run_Seurat_Transfer_PCA.R
│   ├── run_Seurat_Transfer_CCA.R
│   ├── run_Seurat_Transfer_RPCA.R
│   └── run_CIPR.R
├── classic-ML-based/                      # 30+ classic ML tools
│   ├── run_scPred_*.R (19 variants)
│   ├── run_scAnnotatR.R
│   ├── run_CALLR.R
│   ├── run_CHETAH.R
│   ├── run_scClassify.R
│   ├── run_singleCellNet.R
│   ├── run_CellTypist.R
│   ├── run_scAnnotate.R
│   ├── run_scAnno.R
│   ├── run_CaSTLe.R
│   └── run_scID.R
├── DL-based/                              # 4 deep learning tools
│   ├── run_CAMLU.R
│   ├── run_NeuCA.R
│   ├── run_scLearn.R
│   ├── run_scPred_mxnetAdam.R
│   └── NeuCA_helpers/
└── LLM-based/                             # 12 LLM-based scripts
    ├── run_GPTCelltype.R
    ├── run_CASSIA.R
    ├── normalisation_pipeline.R
    └── run_mLLMCelltype_*.R (8 variants + base)
```

Note: `database-based/` and `marker-based/` mirror the same tool set under two configurations
(database-connected vs. matched-reference oracle, see Methods Phase 4 in the manuscript) — each
originally shipped with its own copy of the ~15MB CellMarker2.0 `whole.db` reference database,
which is excluded here. Search for the official CellMarker2.0 distribution to reconstitute it
(the `run_SCSA.R`/`run_scType.R`/etc. scripts in each folder expect a `whole.db` file alongside
them, or a path passed in explicitly).

## Implementation Notes

### Cross-Compatible Seurat Access
Tools should handle both Seurat v4 and v5:
```r
tryCatch({
  data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
}, error = function(e) {
  data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
})
```

### Marker Filtering
**Centralized function** — `prepare_markers()` in `benchmarking_helpers.R` (all tools except SCINA):
```r
filtered <- markers_df %>%
  filter(pct.1 >= 0.10, avg_log2FC > 0) %>%  # NO p_val_adj filter (by design)
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()
```
- `p_val_adj` is **deliberately excluded**: in hard/rare scenarios, low power drives p_val_adj → 1.0 even for genuine markers. Filtering on it would silently empty marker lists.
- Thresholds: `pct.1 ≥ 0.10`, `avg_log2FC > 0`, **top 20 per cell type**.

**Exception — SCINA** (`run_SCINA.R`): Uses its own inline filtering (`avg_log2FC ≥ 0.5`, `p_val_adj < 0.05`, `pct.1 ≥ 0.15`, top 20). Does not call `prepare_markers()`. Results for SCINA use stricter, p_val_adj-filtered markers.

### Fold Caching
`get_fold_data_cached()` stores fold splits and markers in a global `fold_cache` variable to avoid recomputing `FindAllMarkers()` across tool runs. The cache persists across tools within a session but should be cleared (`fold_cache <- list()`) when changing datasets.

### Confidence Scores
- Use tool-specific confidence measures when available (e.g., SingleR delta, Seurat prediction.score.max)
- Default to 1.0 for successful predictions, 0.0 for "Unknown"
- Handle NA values appropriately

## Known Issues

### Critical (affect metric validity)
1. **Garnett ground-truth leakage in fallback**: When `cluster_ext_type` column is absent after `classify_cells()`, Garnett falls back to mapping `garnett_cluster` → cell type via majority vote on `seurat_test$Ground_Truth_Celltype`. This uses test set labels to construct predictions, artificially inflating accuracy. Occurs silently with a cat() message but no warning.
3. **SCINA case mismatch**: SCINA outputs `"unknown"` (lowercase) for unassigned cells. `calculate_metrics()` only filters `"Unknown"` (capitalized), so SCINA's unassigned cells are treated as valid predictions. `unassigned_rate` will always report 0 for SCINA regardless of true unassignment rate.

### High (affect interpretation)
4. **scCATCH cluster vector from ground truth**: When `Ground_Truth_Celltype` is present, scCATCH uses it as the cluster vector for `createscCATCH()`. This means scCATCH "clusters" are the true cell type labels, giving it unfair access to ground truth structure during annotation.
5. **scType network dependency**: Sources three GitHub URLs (`IanevskiAleksandr/sc-type`) on every function call. Fails offline; introduces variability if the remote files change between folds.
6. **scType fallback uses different filtering**: The `run_sctype()` fallback path (direct `sctype_score`) uses hardcoded `0.5/0.05/0.15` thresholds instead of `prepare_markers()`. If the wrapper fails, filtering strategy silently changes.

### Low (annotations/portability)
7. **SCSA conda path**: Falls back to `~/miniforge3/bin/conda` when the resolved `conda` binary looks like a homebrew shim — requires a miniforge3 install at the default location on other machines.
8. **CellAssign TensorFlow memory**: `peakRAM` only captures R-side memory; TensorFlow's memory pool is separate and untracked. Memory stats for CellAssign are systematically underreported.
9. **CellAssign TensorFlow compatibility**: Monkey-patches `tf.reshape()` for compatibility; may still fail with newer TF versions.
10. **CIPR marker recomputation**: Ignores passed markers; computes its own from test data (by design).
11. **Global fold_cache**: Not automatically cleared between datasets — must manually run `fold_cache <- list()`.

## Best Practices

1. **Minimal Implementation**: Keep tool functions focused on core algorithm logic
2. **Error Handling**: Always provide graceful fallbacks to "Unknown" predictions (capital U — the framework's sentinel value)
3. **Memory Tracking**: Wrap core algorithm in `peakRAM::peakRAM()` for memory measurement
4. **Consistency**: Follow established patterns for data access and result formatting
5. **Marker Filtering**: Call `prepare_markers(markers)` from `benchmarking_helpers.R`; do not implement inline filtering
6. **Return format**: `predictions`, `true_labels`, `cell_ids` must all have the same length as `ncol(seurat_test)` — cluster-level tools must expand to cell level before returning
7. **Unknown sentinel**: Use `"Unknown"` (capital U) for unassigned cells — lowercase `"unknown"` is not filtered by `calculate_metrics()`
8. **Testing**: Validate that functions work with the benchmarking framework via `run_single_tool_cv()`
