# Memory Tracking & Marker Usage Audit

**Last updated:** 2026-03-16

## Overview

This document records:
1. Which annotation tools use the `markers` parameter (from `FindAllMarkers`) and why
2. Whether each tool has memory tracking in its main algorithm/tool call section
3. What memory tracking package is used

Memory tracking covers the **core algorithm call only** — Seurat-to-SCE/h5ad data conversions are excluded.

---

## Marker Usage

### Why Markers Are Needed Per Tool Category

**Marker-Based Tools** use the `markers` DataFrame as their primary annotation signal. They have no reference expression data, so they rely entirely on cluster differential expression to identify cell types.

**LLM-Based Tools** pass marker gene lists as text context to the LLM. The LLM reasons semantically over the genes to infer cell type names — the marker list *is* the annotation input.

**Reference-Based Tools (most)** do not need markers because they compare query expression profiles directly against a reference dataset. The `markers` parameter is accepted for interface consistency only.

**Classic ML / DL Tools** train classifiers or embed cells from expression data directly, so markers are not required.

---

### Marker Usage Table

#### Marker-Based Tools (all 10 use markers)

| Tool | How Markers Are Used |
|------|---------------------|
| scCATCH | Converted to scCATCH custom marker database format (species/celltype/gene columns); internal `findmarkergene()` cross-references against it |
| SCINA | Converted to named list of gene signature vectors (one per cell type); used as probabilistic classifiers |
| scType | Converted to Excel-format custom marker DB (positive + negative markers per cell type); scored via `sctype_score()` |
| Garnett | Written to a Garnett-format marker file (text); used to train a cell classifier via `train_cell_classifier()` |
| clustifyr (hyper) | Converted to wide-format ranked gene list; hypergeometric enrichment test against query cluster genes |
| clustifyr (jaccard) | Same as above; Jaccard similarity instead of hypergeometric |
| ScInfeR | Converted to (celltype, marker, weight=1) DataFrame; kNN-based inference from weighted marker scores |
| CellAssign | Top 5 markers per cell type → binary marker gene matrix; probabilistic EM assignment via `cellassign()` |
| scSorter | Converted to (Type, Marker, Weight) annotation DataFrame; semi-supervised cell sorting |
| SCSA | Written to CSV file; passed via command-line to Python SCSA script for database-driven annotation |

#### Reference-Based Tools

| Tool | Uses Markers? | Notes |
|------|---------------|-------|
| SingleR | **YES** | Marker list passed to `SingleR(..., genes = marker_list)` to restrict DE-based scoring to known cell type markers |
| Seurat Transfer (PCA) | No | Interface consistency only; PCA-based transfer learning from expression |
| Seurat Transfer (CCA) | No | Interface consistency only; CCA-based transfer learning |
| Seurat Transfer (RPCA) | No | Interface consistency only; RPCA-based transfer learning |
| scmap_cell | No | Interface consistency only; k-nearest neighbor cell mapping |
| scmap_cluster | No | Interface consistency only; cluster-level projection mapping |
| scibetR | No | Interface consistency only; TPM-based classification via `SciBet_R()` |
| CIPR | **Recomputed internally** | Does NOT use the passed `markers`. Calls `FindAllMarkers(seurat_test, ...)` internally on test data, then compares to training reference expression. Cluster-based design requires test-set markers. |

#### Classic ML-Based Tools

| Tool | Uses Markers? | Notes |
|------|---------------|-------|
| scAnnotatR | **YES** | Markers used to train cell-type-specific classifiers via `train_classifier()` |
| CALLR | **YES** | Converted to a cell_types × top_markers matrix for representative cell selection |
| scPred (all 19 variants) | No | Interface consistency only; train SVM/RF/GLM etc. from PCA embeddings |
| CellTypist | No | Interface consistency only; logistic regression from expression |
| scID | No | Interface consistency only; NMF-based reference matching |
| scAnnotate | No | Interface consistency only |
| scAnno | No | Interface consistency only |
| scClassify | No | Interface consistency only |
| CHETAH | No | Interface consistency only |
| singleCellNet | No | Interface consistency only |
| CaSTLe | No | Interface consistency only; mutual information feature selection from expression |

#### DL-Based Tools

All DL-based tools accept `markers` for interface consistency but do not use them. They train neural networks or embed cells from raw expression data — pre-computed DE markers are neither required nor used.

| Tool | Uses Markers? | Algorithm Type | Why Markers Are Not Needed |
|------|:------------:|----------------|---------------------------|
| CAMLU (R) | No | Deep learning (R wrapper) | Trains from expression profiles |
| NeuCA (R) | No | Deep learning (R wrapper) | Trains from expression profiles |
| scLearn (R) | No | Deep learning (R wrapper) | Trains from expression profiles |
| scPred_mxnetAdam (R) | No | MXNet neural network | Trains from PCA embeddings |
| scDeepSort | No | GNN classifier | Learns graph-structured representations from expression |
| scHash | No | Hash-based DL | Learns hashed representations from expression profiles |
| TOSICA | No | Pathway attention transformer | Uses Reactome GMT gene sets internally; not DE markers |
| scMMT | No | Multi-modal transformer | Learns from raw/normalized expression end-to-end |
| scnym | No | Adversarial domain adaptation | Domain-adapted training on full expression |
| CIForm | No | Consensus integration (transformer) | Feature integration from expression profiles |
| scBalance | No | Dropout NN with class balancing | Weighted NN trained on expression |
| Cell_BLAST | No | DIRECTi VAE + BLAST query | Variational model learns latent space; NN search at inference |
| SCANVI | No | Semi-supervised VAE (scvi-tools) | scVI reference + SCANVI fine-tuning on full expression |
| MARS | No | — | Interface consistency only |
| mtANN | No | — | Interface consistency only |
| scDeepinsight | No | — | Interface consistency only |

#### LLM-Based Tools (all 10 use markers)

| Tool | How Markers Are Used |
|------|---------------------|
| GPTCelltype | Passed directly to `gptcelltype(markers, ...)` |
| CASSIA | Passed to `runCASSIA_pipeline(marker = markers, ...)` |
| mLLMCelltype | Passed to `interactive_consensus_annotation(markers, ...)` |
| mLLMCelltype_claude_sonnet4.5 | Same as above |
| mLLMCelltype_deepseek_v3.1_terminus | Same as above |
| mLLMCelltype_gemini_pro2.5 | Same as above |
| mLLMCelltype_gpt5 | Same as above |
| mLLMCelltype_grok4_fast | Same as above |
| mLLMCelltype_llama4_maverick | Same as above |
| mLLMCelltype_qwen3_max | Same as above |

---

## Memory Tracking Audit

### Packages Used

| Package | Description | Used By |
|---------|-------------|---------|
| `peakRAM::peakRAM()` | R package; wraps the core algorithm call and polls actual process RSS (resident set size) at OS level; extracts `$Peak_RAM_Used_MiB[1]` and `$Elapsed_Time_sec[1]` | All R tools |
| `memory_profiler.memory_usage()` | Python package; polls process RSS at configurable intervals (`interval=0.1`); `include_children=True` captures spawned threads/processes | DL-based Python tools, CellTypist |
| Python stdout parsing | R reads `PEAK_SYSTEM_MEMORY_MB:` and `RUNTIME_SECS:` strings from subprocess stdout | SCSA |

### Standard Patterns

**R tools (`peakRAM::peakRAM()`):**
```r
runtime_secs          <- NA
peak_system_memory_mb <- NA

if (!requireNamespace("peakRAM", quietly = TRUE)) {
  warning("peakRAM package not available for memory/time tracking")
  start_time <- Sys.time()
  result <- tool_function(...)
  runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
} else {
  library(peakRAM)
  peakRAM_result <- peakRAM::peakRAM({
    result <- tool_function(...)
  })
  runtime_secs          <- peakRAM_result$Elapsed_Time_sec[1]
  peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
}
```

The tool returns `runtime_secs` and `peak_system_memory_mb` in its output list. `benchmarking_helpers.R` uses `fold_result$runtime_secs` directly (with outer wall-clock fallback if NULL/NA) and collects `fold_result$peak_system_memory_mb`.

**Python DL tools (`memory_profiler`):**
```python
import time
from memory_profiler import memory_usage

_timing = [0.0]

def _run_method():
    _t0 = time.perf_counter()
    # ... model setup, training, prediction ...
    _timing[0] = time.perf_counter() - _t0

mem_list, _ = memory_usage((_run_method, [], {}), interval=0.1,
                           include_children=True, retval=True)
peak_system_memory_mb = max(mem_list)
method_walltime = _timing[0]

results = {
    'predictions': ...,
    'confidence_scores': ...,
    'peak_system_memory_mb': peak_system_memory_mb,
    'runtime_secs': method_walltime,
    ...
}
```

> **Measurement window:** model setup → training → prediction only. Excludes: conda env startup, Python import time, h5ad file loading, and Seurat→h5ad R-side conversion. `include_children=True` ensures spawned threads (e.g. DataLoader workers) are included in the RSS measurement.

**SCSA (`memory_profiler` at `__main__` level):**
```python
import time
try:
    from memory_profiler import memory_usage
    _has_mp = True
except ImportError:
    _has_mp = False

_start = time.perf_counter()
if _has_mp:
    _mem = memory_usage(lambda: p.run_cmd(args), interval=0.05)
    _peak_mb = max(_mem)
else:
    p.run_cmd(args)
    _peak_mb = float('nan')
_runtime = time.perf_counter() - _start

print(f"PEAK_SYSTEM_MEMORY_MB:{_peak_mb:.2f}")
print(f"RUNTIME_SECS:{_runtime:.4f}")
```

`run_SCSA.R` parses both `PEAK_SYSTEM_MEMORY_MB:` and `RUNTIME_SECS:` from subprocess stdout via `grep()` + `sub()`.

**CellTypist (`memory_profiler` in `celltypist_helper.py`):**

Core annotation (model load, train, predict) is wrapped in `run_core_algorithm()`. `memory_usage()` wraps this function; `time.perf_counter()` provides wall-clock runtime. Both values are emitted to stdout and parsed by `run_CellTypist.R`.

---

### Memory Tracking Per Tool

#### Marker-Based Tools

| Tool | Memory Tracking | Package |
|------|----------------|---------|
| scCATCH | YES | peakRAM |
| SCINA | YES | peakRAM |
| scType | YES | peakRAM (×3: wrapper + 2 fallbacks) |
| Garnett | YES | peakRAM (×2: train + classify) |
| clustifyr_hyper | YES | peakRAM |
| clustifyr_jaccard | YES | peakRAM |
| ScInfeR | YES | peakRAM |
| CellAssign | YES | peakRAM |
| scSorter | YES | peakRAM |
| SCSA | YES | Python stdout (PEAK_SYSTEM_MEMORY_MB: + RUNTIME_SECS:) via memory_profiler at `__main__` level |

#### Reference-Based Tools

| Tool | Memory Tracking | Package |
|------|----------------|---------|
| SingleR | YES | peakRAM |
| Seurat Transfer PCA | YES | peakRAM |
| Seurat Transfer CCA | YES | peakRAM |
| Seurat Transfer RPCA | YES | peakRAM |
| scmap_cell | YES | peakRAM |
| scmap_cluster | YES | peakRAM |
| scibetR | YES | peakRAM |
| CIPR | YES | peakRAM |

#### Classic ML-Based Tools

| Tool | Memory Tracking | Package |
|------|----------------|---------|
| scPred_glm | YES | peakRAM |
| scPred_earth | YES | peakRAM |
| scPred_glmboost | YES | peakRAM |
| scPred_svmPoly | YES | peakRAM |
| scPred_mlp | YES | peakRAM |
| scPred_bayesglm | YES | peakRAM |
| scPred_multinom | YES | peakRAM |
| scPred_avNNet | YES | peakRAM |
| scPred_adaboost | YES | peakRAM |
| scPred_nb | YES | peakRAM |
| scPred_rf | YES | peakRAM |
| scPred_knn | YES | peakRAM |
| scPred_regLogistic | YES | peakRAM |
| scPred_glmnet | YES | peakRAM |
| scPred_nnet | YES | peakRAM |
| scPred_lda | YES | peakRAM |
| scPred_svmLinear | YES | peakRAM |
| scPred_xgbTree | YES | peakRAM |
| scPred_mxnetAdam | YES | peakRAM |
| CaSTLe | YES | peakRAM |
| scClassify | YES | peakRAM |
| scAnnotatR | YES | peakRAM |
| CALLR | YES | peakRAM (×2: preprocess + callr; runtimes summed, peaks max'd) |
| CHETAH | YES | peakRAM |
| singleCellNet | YES | peakRAM |
| scID | YES | peakRAM |
| scAnno | YES | peakRAM |
| scAnnotate | YES | peakRAM |
| CellTypist | YES (Python side) | memory_profiler (celltypist_helper.py) |

#### DL-Based Tools

R wrappers run on the primary machine using `peakRAM`. Python tools run on a separate machine inside isolated conda subprocesses and use `memory_profiler.memory_usage()`. MARS, mtANN, and scDeepinsight have not yet been audited on the remote machine.

| Tool | Memory Tracking | Package | Notes |
|------|----------------|---------|-------|
| CAMLU (R wrapper) | YES | peakRAM | |
| NeuCA (R wrapper) | YES | peakRAM | |
| scLearn (R wrapper) | YES | peakRAM | |
| scDeepSort | YES | memory_profiler | Core algorithm only (model setup → predict) |
| scHash | YES | memory_profiler | Core algorithm only; enforces float32 |
| TOSICA | YES | memory_profiler | Core algorithm only |
| scMMT | YES | memory_profiler | Core algorithm only |
| scnym | YES | memory_profiler | Core algorithm only |
| CIForm | YES | memory_profiler | Core algorithm only |
| scBalance | YES | memory_profiler | Core algorithm only |
| Cell_BLAST | YES | memory_profiler | Core algorithm only |
| SCANVI | YES | memory_profiler | Core algorithm only |
| MARS | Not audited | — | Remote machine, not yet verified |
| mtANN | Not audited | — | Remote machine, not yet verified |
| scDeepinsight | Not audited | — | Remote machine, not yet verified |

#### LLM-Based Tools

| Tool | Memory Tracking | Package |
|------|----------------|---------|
| GPTCelltype | YES | peakRAM |
| CASSIA | YES | peakRAM |
| mLLMCelltype | YES | peakRAM |
| mLLMCelltype_claude_sonnet4.5 | YES | peakRAM |
| mLLMCelltype_deepseek_v3.1_terminus | YES | peakRAM |
| mLLMCelltype_gemini_pro2.5 | YES | peakRAM |
| mLLMCelltype_gpt5 | YES | peakRAM |
| mLLMCelltype_grok4_fast | YES | peakRAM |
| mLLMCelltype_llama4_maverick | YES | peakRAM |
| mLLMCelltype_qwen3_max | YES | peakRAM |

---

## DL-Based Python Tools — Implementation Notes

These tools run on a separate machine in isolated conda subprocesses. `memory_profiler.memory_usage()` wraps only the core algorithm (`_run_method()`), covering model setup → training → prediction. This deliberately **excludes** conda environment startup, Python import time, h5ad file loading, and Seurat→h5ad preprocessing. `include_children=True` ensures memory from spawned threads (e.g. DataLoader workers) is captured. Runtime is measured with `time.perf_counter()` inside `_run_method()`.

| Tool | Covers h5ad Load | Covers Training | Covers Prediction |
|------|:----------------:|:---------------:|:-----------------:|
| scDeepSort | No | Yes | Yes |
| scHash | No | Yes | Yes |
| TOSICA | No | Yes | Yes |
| scMMT | No | Yes | Yes |
| scnym | No | Yes | Yes |
| CIForm | No | Yes | Yes |
| scBalance | No | Yes | Yes |
| Cell_BLAST | No | Yes | Yes |
| SCANVI | No | Yes | Yes |

**Per-tool notes:**

- **scDeepSort** — Uses CSV format (not h5ad) for gene/cell data; `num_neighbors=100, n_layers=1` to avoid over-smoothing; GPU disabled (`gpu_id=-1`)
- **scHash** — Enforces float32 (`adata.X = adata.X.astype(np.float32)`) due to PyTorch Linear layer requirements
- **TOSICA** — Uses Reactome GMT gene sets for pathway definitions; `laten=True` enables latent path bypass; returns per-cell probability scores
- **scMMT** — Complex preprocessing (normalize_total → log1p → scale); aggressive NaN/Inf cleaning and value clipping; `gene_normalize=True` required
- **scnym** — Requires CPM normalization (1e6, not 1e4); uses `scNym_confidence` column for native confidence scores; trains with explicit `train`/`test` domain labels
- **CIForm** — Uses `run_CIForm_helper.py`; case-insensitive label normalization (lowercase); fixed seed s=1024, 50 training epochs
- **scBalance** — Enforces float32 for PyTorch compatibility; weighted sampling for class imbalance; DataFrame interface (cells × genes)
- **Cell_BLAST** — Requires raw count data (reverses log-normalization via `expm1`); variable gene selection via DIRECTi; majority voting at inference (`min_hits=1, majority_threshold=0.5`)
- **SCANVI** — Two-stage pipeline: SCVI reference model → SCANVI fine-tuning; online query mapping with online update training; confidence from max class probability per cell
