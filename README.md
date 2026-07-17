# scRNA-seq Cell-Type Annotation Benchmark — Code Archive

This is a code-only reference archive extracted from the working repository behind a thesis
benchmarking ~58 single-cell RNA-seq cell-type annotation tools across five experimental phases
(controlled synthetic simulation through to real-world, deployment-realistic conditions). It
exists so readers of the manuscript can see how the analysis was actually implemented.

**This repository contains no data.** Raw and simulated single-cell datasets, preprocessed
Seurat/AnnData objects, the CellMarker2.0 reference database, benchmark result tables, figures,
and Python virtual environments are all excluded — see [What's excluded](#whats-excluded) below.
Only the R/Python source code and analysis scripts are archived.

## Pipeline overview

The study proceeds in three stages, mirrored by the three top-level folders here:

1. **[`01-data-generation/`](01-data-generation/)** — generates the Phase 1 synthetic (Taguchi
   L9(3⁴) orthogonal design) data and preprocesses the real and cross-platform datasets used in
   later phases, including the profiling/Borda-consensus procedure used to match real datasets to
   L9 scenarios.
2. **[`02-benchmarking-pipeline/`](02-benchmarking-pipeline/)** — the shared benchmarking
   framework and all ~58 tool implementations (marker-based, reference-based, classic ML, deep
   learning, and LLM-based), used across every phase.
3. **[`03-results-analysis/`](03-results-analysis/)** — the statistics/figure-generation scripts
   for each reported phase.

## Manuscript phase mapping

The manuscript narrates five phases; this codebase's `03-results-analysis/` folders are named
descriptively rather than by their original internal phase codes. Some folders support more than
one phase's discussion (marked *supplement* below), since the underlying tool comparisons are
shared across phases.

| Manuscript phase | Description | `03-results-analysis/` folder(s) |
|---|---|---|
| **Phase 1** | Factorial Splatter simulation, Taguchi L9(3⁴) orthogonal design, oracle inputs — main-effect attribution | `synthetic-L9-design/`, plus `marker-synthetic-supplement/` and `foundation-models-synthetic/` |
| **Phase 2** | Nine within-platform real datasets — tests whether Phase 1 structure-to-accuracy relationships survive on real biology | `real-data-validation/`, plus `marker-real-supplement/` and `foundation-models-real/` |
| **Phase 3** | Cross-platform transfer across three technology blocks (CellBench, Pancreas, PBMCbench) | `cross-platform-transfer/` |
| **Phase 4** | Database-connected marker tools and LLM tools, ontology-aware 0/0.5/1 scoring panel | `database-marker-tools/`, `llm-ontology-scoring/` |
| **Phase 5** | Foundation models (transformer-based); real-pretrained arm is the reported result | `foundation-models-real/` (reported), `foundation-models-synthetic/` (supplementary) |

Each folder contains that phase's `Plots_Analysis.R` (statistics + figure generation), sourcing
the shared `03-results-analysis/cluster_dependent_helper.R`.

The nine real datasets used in Phase 2+ (selected from a larger candidate pool via a
Borda-consensus procedure over three dataset-profiling schemes — see
`01-data-generation/dataset-profiling/`) are: Darmanis-Brain-2015, Marques-Brain-2016,
Nowakowski-Cortex-2017, Grün-Pancreas-2016, Tabula Muris-FACS-3k, He-Skin-2020,
Zhao-Immune-Fine-2020, Zheng-ZhengSort-5cl-2017, and MacParland-Liver-Broad-2018.

## What's excluded

- **Raw/processed data**: `.rds`, `.h5ad`, `.mtx`, `.tsv`/`.csv` data tables, `.RData`
- **Figures/plots**: all `.png`/`.pdf` outputs (the scripts that generate them are included)
- **Logs**: run logs, LLM API transcripts
- **CellMarker2.0 reference database** (`whole.db`, ~15MB) — used by the marker-based/
  database-based tool categories; search for the official CellMarker2.0 distribution to
  reconstitute it
- **Python virtual environments** — see `02-benchmarking-pipeline/docs/` for setup notes
- **Superseded pilot work**: early trial phases and V1 dataset-preparation passes that predate
  the design reported in the manuscript are not included
- **Manuscript text, drafts, and build tooling** — this is a code-only archive; see the published
  manuscript for the narrative

## Reproducing a run

This archive is a reference, not a turnkey pipeline — it assumes you'll supply your own copies of
the source datasets and the CellMarker2.0 database, and reconstruct the R/Python environments
described in `02-benchmarking-pipeline/docs/`. Scripts that read data use relative paths or
`Sys.getenv()`-configurable paths with sensible defaults (documented inline where non-obvious);
point them at your own data directory structure following the phase layout above.

## Structure

```
├── 01-data-generation/         # synthetic simulation, real/cross-platform preprocessing, loaders, profiling
├── 02-benchmarking-pipeline/   # shared framework + all tool implementations
└── 03-results-analysis/        # per-phase statistics and figure generation
```

See each subfolder's own README for details.
