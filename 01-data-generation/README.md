# Data Generation

Scripts that generate the Phase 1 synthetic data and preprocess the real/cross-platform datasets
used in later phases. No data is included — these are the generation/preprocessing scripts only.

## `synthetic-simulation/`

Generates the Taguchi L9(3⁴) orthogonal-array synthetic benchmark (Phase 1) using Splatter.
`L9_Generator.R` produces 9 scenarios × 3 replicates; `L9_Sanity_Check.R` validates the output.

## `real-data-preprocessing/`

One `*_Preprocessing.R` (or `.ipynb`) script per real candidate dataset — QC, normalization, and
standardization into the `Ground_Truth_Celltype`-labeled format the benchmarking framework
expects. This is the full candidate pool profiled and Borda-matched to L9 scenarios (see
`dataset-profiling/`); only nine of these datasets were ultimately selected for the reported Phase
2 panel (see the top-level README's phase-mapping table).

## `cross-platform-preprocessing/`

Preprocessing scripts for the Phase 3 cross-platform transfer datasets (CellBench across three
sequencing protocols, PBMCBench across five, Baron-Pancreas, Segerstolpe-Pancreas).

## `data-loaders/`

Converts the raw simulation/preprocessing output into Seurat (`.rds`) and AnnData (`.h5ad`)
objects so both R- and Python-based tools receive comparable inputs. Filenames indicate which
stage they load for (`_P4_L9`, `_P5`, `_P9`).

## `dataset-profiling/`

Computes the dataset-profile metrics (cell count, class balance, DE strength, kNN purity, etc.)
used to match candidate real datasets to their nearest L9 scenario via Borda-consensus over three
complementary profiling schemes (10-metric profile, 4-axis PCA, real-data PCA). Produces the
`Profiling_L9.csv` / `Profiling_Real.csv` inputs consumed by several `03-results-analysis/`
scripts (not included here — regenerate by running these notebooks against your own data).
