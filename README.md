# scRNA-seq Cell-Type Annotation Benchmark — Code Archive

This is a code-only reference archive extracted from the working repository behind a thesis
benchmarking ~63 single-cell RNA-seq cell-type annotation tools across five experimental phases
(controlled synthetic simulation through to real-world, deployment-realistic conditions). 

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
descriptively rather than by their original internal phase codes. 

| Manuscript phase | Description | `03-results-analysis/` folder(s) |
|---|---|---|
| **Phase 1** | Factorial Splatter simulation, Taguchi L9(3⁴) orthogonal design, oracle inputs — main-effect attribution | `synthetic-L9-design/`|
| **Phase 2** | Nine within-platform real datasets — tests whether Phase 1 structure-to-accuracy relationships survive on real biology | `real-data-validation/`|
| **Phase 3** | Cross-platform transfer across three technology blocks (CellBench, Pancreas, PBMCbench) | `cross-platform-transfer/` |
| **Phase 4** | Database-connected marker tools and LLM tools, ontology-aware 0/0.5/1 scoring panel | `database-marker-tools/`, `llm-ontology-scoring/` |
| **Phase 5** | Foundation models (transformer-based); real-pretrained arm is the reported result | `foundation-models-real/` |

The nine real datasets used in Phase 2+ (selected from a larger candidate pool via a
Borda-consensus procedure over three dataset-profiling schemes — see
`01-data-generation/dataset-profiling/`) are: Darmanis-Brain-2015, Marques-Brain-2016,
Nowakowski-Cortex-2017, Grün-Pancreas-2016, Tabula Muris-FACS-3k, He-Skin-2020,
Zhao-Immune-Fine-2020, Zheng-ZhengSort-5cl-2017, and MacParland-Liver-Broad-2018.
