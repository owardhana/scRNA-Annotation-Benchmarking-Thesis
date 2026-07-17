# Marker Filtering Standardization

**Last updated:** 2026-03-25

---

## Overview

All marker-using annotation tools receive markers pre-filtered by a single centralized function, `prepare_markers()`, defined in `benchmarking_helpers.R`. This replaces per-tool inline filtering that previously caused inconsistencies.

---

## Standard Filtering Parameters

Implemented via `prepare_markers(markers_df, top_n = 20)` in `benchmarking_helpers.R`:

```r
prepare_markers <- function(markers_df, top_n = 20) {
  filtered <- markers_df %>%
    filter(
      pct.1 >= 0.10,       # Gene expressed in ≥10% of target cells
      avg_log2FC > 0        # Positive markers only
    ) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()
  return(filtered)
}
```

### Rationale for parameter choices

| Parameter | Value | Reason |
|-----------|-------|--------|
| `avg_log2FC > 0` | positive only | Selects the strongest available signal. No minimum FC threshold — in hard DE (FC ~1.5×) and rare-type scenarios, a 0.5 threshold silently produces empty marker lists. |
| `pct.1 >= 0.10` | ≥10% cells | Removes near-zero-expression genes while tolerating rare cell types with low marker coverage. |
| `p_val_adj` | **excluded** | Not a discovery analysis. With hard DE scenarios and rare types (1–3 training cells), p_val_adj ≈ 1.0 due to low power — not poor marker quality. Filtering by it would eliminate markers for the exact scenarios we need to benchmark. |
| `top_n = 20` | 20 per type | Balances signal richness with computational overhead. Uniform across tools for fair comparison. |

---

## Tools Using `prepare_markers()`

### Marker-Based (10 tools)

| Tool | How markers are used | Notes |
|------|---------------------|-------|
| scCATCH | Custom marker database (species/celltype/gene) | Standard |
| scType | Excel-format marker DB (positive/negative per type) | Standard (main + fallback) |
| Garnett | Garnett-format marker file | Standard |
| SCINA | Named list of gene signature vectors per cell type | Standard; inline `head(..., 20)` hardcoded in `convert_markers_to_scina()` to match |
| clustifyr_hyper | Wide-format ranked gene list; hypergeometric enrichment | Standard |
| clustifyr_jaccard | Same as hyper; Jaccard similarity metric | Standard |
| ScInfeR | (celltype, marker, weight=1) DataFrame | Standard |
| scSorter | (Type, Marker, Weight) annotation DataFrame | Standard |
| CellAssign | Binary marker gene matrix for EM assignment | Standard |
| SCSA | CSV file passed to Python; `-f 0` in CLI (no additional FC threshold) | See exception below |

### Reference-Based (1 tool)

| Tool | How markers are used |
|------|---------------------|
| SingleR | Gene list passed to `SingleR(..., genes = marker_list)` to restrict DE-based scoring |

### Classic ML-Based (2 tools)

| Tool | How markers are used | Notes |
|------|---------------------|-------|
| scAnnotatR | Trains cell-type-specific SVM classifiers via `train_classifier()` | `head(ct_df$gene, 20)` enforced inline after sorting by FC |
| CALLR | Cell_types × top_markers matrix for representative cell selection | `prepare_markers(markers, top_n = 20)` explicit override |

### LLM-Based (10 tools)

All LLM tools call `prepare_markers(markers)` at function entry (inheriting `top_n = 20`), then pass `top_gene_count = 20` to the LLM annotation function as a secondary guardrail.

| Tool | LLM function |
|------|-------------|
| GPTCelltype | `gptcelltype()` |
| CASSIA | `runCASSIA_pipeline()` |
| mLLMCelltype | `interactive_consensus_annotation()` |
| mLLMCelltype_claude_sonnet4.5 | `annotate_cell_types()` |
| mLLMCelltype_deepseek_v3.1_terminus | `annotate_cell_types()` |
| mLLMCelltype_gemini_pro2.5 | `annotate_cell_types()` |
| mLLMCelltype_gpt5 | `annotate_cell_types()` |
| mLLMCelltype_grok4_fast | `annotate_cell_types()` |
| mLLMCelltype_llama4_maverick | `annotate_cell_types()` |
| mLLMCelltype_qwen3_max | `annotate_cell_types()` |

---

## Exception: SCSA

**File**: `run_SCSA.R`

SCSA uses `prepare_markers()` in R (so the standard top-20 / pct.1 / FC > 0 filter applies before writing the CSV), but the Python SCSA command uses `-f 0` (no additional fold change threshold on the Python side). This prevents double-filtering.

SCSA has no `top-N` selection capability in Python — it uses all rows in the input CSV. The top-20 limit is therefore entirely enforced by `prepare_markers()` in R.

```r
markers_to_write <- prepare_markers(markers)   # top 20, pct.1 >= 0.10, FC > 0
# Python CLI:
scsa_args <- c(..., "-f", "0", ...)            # no additional FC filter in Python
```

---

## Complete Filtering Summary Table

| Tool | FC threshold | pct.1 | p_val_adj | top-N | Method |
|------|-------------|-------|-----------|-------|--------|
| scCATCH | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| scType | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` (main + fallback) |
| Garnett | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| SCINA | > 0 | ≥ 0.10 | excluded | 20 | Inline `head(..., 20)` after 0.5/0.05/0.15 filter |
| clustifyr_hyper | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| clustifyr_jaccard | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| ScInfeR | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| scSorter | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| CellAssign | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| SCSA | > 0 (R) / none (Python) | ≥ 0.10 (R only) | excluded | 20 (R only) | `prepare_markers()` + `-f 0` |
| SingleR | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| scAnnotatR | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` + inline `head(..., 20)` |
| CALLR | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers(markers, top_n = 20)` |
| GPTCelltype | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| CASSIA | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` |
| mLLMCelltype (×9) | > 0 | ≥ 0.10 | excluded | 20 | `prepare_markers()` + `top_gene_count = 20` |

---

## Tools That Do NOT Use Markers

These tools accept `markers` for interface consistency but ignore it:

- **Reference-based**: Seurat Transfer (PCA/CCA/RPCA), scmap_cell, scmap_cluster, scibetR
- **CIPR**: recomputes markers internally from test data
- **Classic ML**: scPred (19 variants), CellTypist, scAnnotate, scAnno, scClassify, CHETAH, singleCellNet, CaSTLe, scID
- **DL-based**: all tools (train from expression data)

---

## History

| Date | Change |
|------|--------|
| 2025-01-06 | Initial per-tool standardization (FC ≥ 0.5, p_val_adj < 0.05, pct.1 ≥ 0.15, top 50) |
| 2025-02-13 | Reinstated after brief unfiltered experiment caused tool failures |
| 2026-03-25 | Replaced all per-tool inline filtering with centralized `prepare_markers()`. Removed p_val_adj filter. Changed thresholds to FC > 0 / pct.1 ≥ 0.10. Extended coverage to SingleR, scAnnotatR, CALLR, all LLM tools. Standardized top-N to 20. |
