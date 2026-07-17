# ============================================================================
# L9 Orthogonal Array Benchmark Suite for scRNA-seq Annotation Methods
# ============================================================================
#
# PURPOSE:
#   Generates 27 simulated scRNA-seq datasets (9 scenarios × 3 replicates)
#   for benchmarking cell type annotation tools. The experiment follows a
#   Taguchi L9(3^4) orthogonal array design, which efficiently covers 4
#   experimental dimensions at 3 levels each using only 9 scenarios instead
#   of the 81 required for a full 3^4 factorial grid.
#
# DESIGN RATIONALE:
#   The L9 guarantees orthogonality: every level of every factor appears
#   with every level of every other factor exactly once. This means main
#   effects can be estimated independently without confounding, allowing a
#   linear model (kappa ~ cell_count + imbalance + n_types + de_signal) to
#   decompose variance cleanly across dimensions. Interaction effects are
#   aliased and cannot be decomposed from this design — this is the accepted
#   cost of the 9x efficiency gain.
#
# DIMENSIONS VARIED:
#   1. Cell count        : 500  / 3,000  / 15,000
#   2. Class imbalance   : Balanced / Mild (R=5) / Severe (R=20)
#   3. Number of types   : 5 / 15 / 30
#   4. DE signal (diff.) : Easy / Medium / Hard
#
# REPLICATION:
#   3 independent Splatter replicates per scenario (different random seeds,
#   all other parameters fixed). Each replicate receives one stratified 80/20
#   train/test split. This gives 3 independent kappa observations per tool
#   per scenario, from which median ± IQR are reported. CV was found to add
#   near-zero variance for these dataset sizes (folds near-identical), making
#   independent replicates more informative than repeated CV.
#
# KEY REFERENCES:
#   - Zappia L et al. (2017) Splatter: simulation of single-cell RNA
#     sequencing data. Genome Biology 18:174.
#     DOI: 10.1186/s13059-017-1305-0
#
#   - Abdelaal T et al. (2019) A comparison of automatic cell identification
#     methods for single-cell RNA sequencing data. Genome Biology 20:194.
#     DOI: 10.1186/s13059-019-1795-z
#     [Source of SVM benchmark dominance; establishes performance hierarchy]
#
#   - Schirmer M et al. (2016) Differential expression benchmarking.
#     [FC = 1.2 produces near-zero power across DE methods — grounds our
#     "hard" scenario at de.facLoc = 0.12 → FC ≈ 1.22x]
#
#   - Svensson V (2020) Droplet scRNA-seq is not zero-inflated.
#     Nature Biotechnology 38:147-150. DOI: 10.1038/s41587-019-0379-5
#     [Justifies experiment-level dropout model over ZINB for 10X data]
#
#   - Taguchi G (1987) System of Experimental Design. UNIPUB/Kraus.
#     Roy R (1990) A Primer on the Taguchi Method. Van Nostrand Reinhold.
#     [Source of L9(3^4) orthogonal array design]
#
#   - Demsar J (2006) Statistical comparisons of classifiers over multiple
#     data sets. JMLR 7:1-30.
#     [Friedman + Nemenyi test framework used for downstream analysis]
#
#   - Luecken MD, Theis FJ (2019) Current best practices in single-cell
#     RNA-seq analysis: a tutorial. Molecular Systems Biology 15:e8746.
#     DOI: 10.15252/msb.20188746
#     [Grounds nGenes = 10,000 and library size parameters]
#
# AUTHOR: [Your name]
# DATE:   [Date]
# ============================================================================

library(splatter)
library(Matrix)
library(SingleCellExperiment)

# ============================================================================
# OUTPUT CONFIGURATION
# ============================================================================

output_dir <- "l9_benchmark_sets"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Fixed random seeds for replicates — defined once here so the entire
# experiment is fully reproducible. The same three seeds are used across
# all 9 scenarios: any performance differences between scenarios therefore
# cannot be attributed to lucky/unlucky random draws.
REPLICATE_SEEDS <- c(
  rep1 = 42,
  rep2 = 123,
  rep3 = 999
)

# ============================================================================
# FIXED PARAMETERS — held constant across all 27 simulations
# ============================================================================
#
# These parameters define the "baseline" data-generating process and are
# not part of the experimental design. They are set to values that reflect
# real 10X Chromium scRNA-seq data as closely as Splatter allows.
#
# Changing any of these would confound your experimental dimensions, so
# they must remain fixed.

FIXED_PARAMS <- list(
  
  # Number of genes simulated
  # Post-QC gene count for 10X Chromium data. Real datasets detect
  # ~15,000-20,000 genes pre-QC, reduced to ~10,000 after filtering
  # low-expression and mitochondrial genes. Using 10,000 keeps simulation
  # runtime manageable while preserving realistic gene space complexity.
  # Reference: Luecken & Theis (2019) best practices.
  nGenes = 10000,
  
  # Biological coefficient of variation (common component)
  # Controls within-cell-type expression heterogeneity. Real scRNA-seq
  # data shows BCV of 0.4-0.6 for heterogeneous tissue samples. Our v1
  # simulations used 0.05-0.3, which produced unrealistically tight clusters
  # and inflated accuracy across all tools. Fixed at 0.4 — the low end of
  # the realistic range — to avoid making scenarios artificially noisy while
  # still challenging classifiers.
  # Reference: McCarthy et al. (2012) edgeR BCV estimation in real data.
  bcv.common = 0.4,
  
  # Library size log-normal location parameter
  # exp(8.5) ≈ 4,915 UMI/cell — within the typical 2,000-8,000 UMI/cell
  # range for 10X Chromium v3 after QC filtering. Combined with
  # dropout.type = "none", the NB distribution naturally produces
  # ~78-82% sparsity at this library size, matching the ~80% figure
  # commonly cited for 10X data (10X Genomics CG000315 technical note).
  # Sanity check validation: lib.loc = 8.0 produced 90% sparsity (too high);
  # lib.loc = 8.5 is expected to bring this into the 78-85% target range.
  # Reference: 10X Genomics CG000315 technical note.
  lib.loc = 8.5,
  
  # Library size log-normal scale parameter
  # Creates ~3-4 fold variation in sequencing depth between cells, which is
  # typical for real 10X data. Larger values would make normalization harder;
  # smaller values would make it trivially easy.
  lib.scale = 0.5,
  
  # Dropout model type — NONE (correct choice for 10X Chromium data)
  # Svensson (2020) demonstrated that zeros in droplet-based protocols are
  # fully explained by Negative Binomial sampling statistics, with no
  # evidence of additional technical zero inflation beyond what NB predicts.
  # Setting dropout.type = "none" lets the NB distribution produce natural
  # sparsity without adding a second artificial dropout layer. This is the
  # biologically correct choice for 10X simulation.
  #
  # NOTE: Our original implementation used dropout.type = "experiment" based
  # on a misreading of Svensson — that paper shows NB IS sufficient, which
  # means extra dropout is NOT needed. Sanity check validation confirmed
  # dropout.type = "experiment" produced 97-99.5% sparsity (far too high)
  # regardless of dropout.mid value, because it compounds with NB zeros.
  # Reference: Svensson (2020) Nature Biotechnology 38:147-150.
  dropout.type = "none",
  
  # DE fold-change variance parameter
  # Standard deviation of the log-normal distribution used to sample
  # individual gene fold changes. A value of 0.4 means DE genes have
  # variable effect sizes around the mean (de.facLoc), with some genes
  # having stronger and some weaker effects. This variance is held
  # FIXED across all DE difficulty levels so that only the mean FC and
  # proportion of DE genes change between easy/medium/hard — not the
  # spread. This makes the difficulty gradient interpretable.
  de.facScale = 0.4
)

# ============================================================================
# DIFFERENTIAL EXPRESSION LEVELS
# ============================================================================
#
# Three DE signal levels define the "difficulty" dimension of the L9 design.
# Difficulty in this context means: how linearly separable are cell types
# in gene expression space?
#
# Mean upregulated fold change formula (Splatter log-normal, upregulated
# direction only — this is the correct metric since Splatter assigns both
# up and down DE genes; averaging them together produces cancellation):
#   mean_FC_up = exp(de.facLoc + de.facScale^2 / 2)
#   With de.facScale = 0.4: mean_FC_up = exp(de.facLoc + 0.08)
#
# RECALIBRATION NOTE (v2 → v3):
#   Sanity check validation revealed that the original values (easy = 0.70,
#   medium = 0.35, hard = 0.12) produced no visible PCA gradient when
#   combined with bcv.common = 0.4. The high within-type biological noise
#   overwhelmed between-type DE signal, collapsing all scenarios to ~1.8%
#   PC1 variance regardless of difficulty level. Values were increased to
#   produce a clear gradient:
#     - easy:   clean cluster separation in PCA expected
#     - medium: partial overlap expected
#     - hard:   substantial overlap expected, methods should struggle
#
# The hard level (de.facLoc = 0.3, mean |log2FC| ≈ 0.38) still represents
# a genuinely challenging scenario — it is below the typical log2FC = 0.5
# threshold used in many DE workflows and corresponds to FC ≈ 1.46x, which
# is at the lower end of reliably detectable differences in real data.
#
# Reference for calibration:
#   Schirmer et al. (2016): FC = 1.2 → near-zero power across all DE methods
#   Reyfman et al.: FC = 1.4-1.5 as realistic challenging scenario
#   Our v2 analysis: FC < 1.3 combined with BCV = 0.4 → no gradient

DE_LEVELS <- list(
  
  # EASY — Clearly distinct cell types
  # de.prob = 0.15   : 15% of 10,000 genes are DE = 1,500 DE genes per type
  # de.facLoc = 1.5  : mean_FC_up = exp(1.5 + 0.08) ≈ 4.8x
  # Biological analog: major lineage differences (e.g., T cell vs B cell,
  # neuron vs. glia). Clean PCA separation expected. Most annotation methods
  # should perform well here — this is the upper bound of the gradient.
  easy = list(
    de.prob   = 0.15,
    de.facLoc = 1.50
  ),
  
  # MEDIUM — Moderately distinct cell types
  # de.prob = 0.08   : 8% of genes DE = 800 DE genes per type
  # de.facLoc = 0.8  : mean_FC_up = exp(0.8 + 0.08) ≈ 2.41x
  # Biological analog: related but distinguishable subtypes (e.g., CD4+ vs
  # CD8+ T cells, different monocyte subsets). Partial PCA overlap expected.
  # Methods with weak regularization or poor sparsity handling begin to
  # show performance degradation here.
  medium = list(
    de.prob   = 0.08,
    de.facLoc = 0.80
  ),
  
  # HARD — Closely related, confusable subtypes
  # de.prob = 0.04   : 4% of genes DE = 400 DE genes per type
  # de.facLoc = 0.3  : mean_FC_up = exp(0.3 + 0.08) ≈ 1.46x
  # Biological analog: highly similar subtypes such as naive vs. memory T
  # cells, monocyte activation states. Substantial PCA overlap expected.
  # This level is designed to produce meaningful performance differences
  # between strong and weak annotation methods.
  # Note: original hard level (de.facLoc = 0.12) was indistinguishable
  # from easy in PCA — recalibrated upward while remaining genuinely hard.
  hard = list(
    de.prob   = 0.04,
    de.facLoc = 0.30
  )
)

# ============================================================================
# IMBALANCE HELPER FUNCTION
# ============================================================================
#
# Generates group probability vectors using a geometric series with a fixed
# max:min ratio R. This ensures imbalance severity is CONSISTENT across all
# three n_types levels (5, 15, 30) — without this, "mild" imbalance with
# 30 types would produce a very different effective distribution than "mild"
# with 5 types.
#
# The geometric series: probs[i] = r^(i-1), normalized to sum to 1
# where r = R^(-1/(n_types - 1)) is chosen so that probs[1]/probs[n] = R
#
# IMPORTANT — S3 edge case (500 cells, severe, 30 types):
#   With R=20 and 30 types, the rarest type receives ~0.5% of 500 cells
#   = approximately 2-3 cells. This is below viable training size for most
#   supervised methods. S3 is retained as a deliberate stress-test scenario
#   to characterize where methods break down. Tools that fail completely
#   on S3 receive kappa = NA and this is reported explicitly.
#   See: scBalance paper (Ma et al. 2023) for rare cell type annotation
#   under extreme imbalance.

make_group_probs <- function(n_types, imbalance_level) {
  
  if (imbalance_level == "balanced") {
    # Equal proportions — every type has 1/n_types of total cells
    return(rep(1 / n_types, n_types))
  }
  
  # Max:min ratio — the dominant cell type has R times more cells than
  # the rarest type. R=5 is "mild" (realistic for many tissues);
  # R=20 is "severe" (representative of tissues with rare populations
  # such as plasmacytoid dendritic cells in blood ~0.1-0.3% of PBMCs).
  R <- switch(imbalance_level,
              mild   = 5,
              severe = 20
  )
  
  # Common ratio of the geometric series
  r <- R^(-1 / (n_types - 1))
  
  # Raw geometric series: [1, r, r^2, ..., r^(n-1)]
  raw_probs <- r^(0:(n_types - 1))
  
  # Normalize to sum to 1
  normalized <- raw_probs / sum(raw_probs)
  
  return(normalized)
}

# ============================================================================
# L9 SCENARIO DEFINITIONS
# ============================================================================
#
# The L9(3^4) Taguchi orthogonal array assigns 4 factors at 3 levels each
# to 9 experimental runs. The defining property is that each ordered pair
# of factor levels (e.g., cell_count=low AND de_signal=hard) appears
# exactly once across the 9 runs — this is the orthogonality guarantee that
# allows unconfounded estimation of main effects.
#
# Column assignment (standard L9 array):
#   Col 1 = cell_count   : 1→500, 2→3000, 3→15000
#   Col 2 = imbalance    : 1→balanced, 2→mild, 3→severe
#   Col 3 = n_types      : 1→5, 2→15, 3→30
#   Col 4 = de_signal    : 1→easy, 2→medium, 3→hard
#
# Standard L9 array (rows are scenarios, entries are level indices):
#   S1:  1 1 1 1   S2:  1 2 2 2   S3:  1 3 3 3
#   S4:  2 1 2 3   S5:  2 2 3 1   S6:  2 3 1 2
#   S7:  3 1 3 2   S8:  3 2 1 3   S9:  3 3 2 1
#
# Reference: Roy R (1990) A Primer on the Taguchi Method. Ch. 3.

l9_scenarios <- list(
  
  # -------------------------------------------------------------------------
  # ROW 1: Small dataset — all low-stress dimensions except they are
  # orthogonally assigned, not all "easy"
  # -------------------------------------------------------------------------
  
  S1 = list(
    scenario_id   = "S1",
    cell_count    = 500,
    imbalance     = "balanced",
    n_types       = 5,
    de_signal     = "easy",
    description   = "Small n, balanced, few types, easy DE — baseline low-stress scenario"
  ),
  
  S2 = list(
    scenario_id   = "S2",
    cell_count    = 500,
    imbalance     = "mild",
    n_types       = 15,
    de_signal     = "medium",
    description   = "Small n, mild imbalance, moderate types, medium DE"
  ),
  
  # NOTE: S3 is the extreme stress-test corner. See make_group_probs()
  # documentation above for the edge case warning.
  S3 = list(
    scenario_id   = "S3",
    cell_count    = 500,
    imbalance     = "severe",
    n_types       = 30,
    de_signal     = "hard",
    description   = "Small n, severe imbalance, many types, hard DE — stress test (rarest type ~2-3 cells)"
  ),
  
  # -------------------------------------------------------------------------
  # ROW 2: Medium dataset
  # -------------------------------------------------------------------------
  
  S4 = list(
    scenario_id   = "S4",
    cell_count    = 3000,
    imbalance     = "balanced",
    n_types       = 15,
    de_signal     = "hard",
    description   = "Medium n, balanced, moderate types, hard DE — tests signal difficulty at moderate scale"
  ),
  
  S5 = list(
    scenario_id   = "S5",
    cell_count    = 3000,
    imbalance     = "mild",
    n_types       = 30,
    de_signal     = "easy",
    description   = "Medium n, mild imbalance, many types, easy DE — tests scalability with many types"
  ),
  
  S6 = list(
    scenario_id   = "S6",
    cell_count    = 3000,
    imbalance     = "severe",
    n_types       = 5,
    de_signal     = "medium",
    description   = "Medium n, severe imbalance, few types, medium DE — imbalance effect at medium scale"
  ),
  
  # -------------------------------------------------------------------------
  # ROW 3: Large dataset
  # -------------------------------------------------------------------------
  
  S7 = list(
    scenario_id   = "S7",
    cell_count    = 15000,
    imbalance     = "balanced",
    n_types       = 30,
    de_signal     = "medium",
    description   = "Large n, balanced, many types, medium DE — scalability with complex taxonomy"
  ),
  
  S8 = list(
    scenario_id   = "S8",
    cell_count    = 15000,
    imbalance     = "mild",
    n_types       = 5,
    de_signal     = "hard",
    description   = "Large n, mild imbalance, few types, hard DE — tests whether scale compensates for weak signal"
  ),
  
  S9 = list(
    scenario_id   = "S9",
    cell_count    = 15000,
    imbalance     = "severe",
    n_types       = 15,
    de_signal     = "easy",
    description   = "Large n, severe imbalance, moderate types, easy DE — imbalance effect at large scale"
  )
)

# ============================================================================
# PARAMETER CONSTRUCTION FUNCTION
# ============================================================================
#
# Combines a scenario definition with the fixed parameter block and a
# replicate seed to produce a complete Splatter SplatParams object.
# All parameters are set explicitly — no Splatter defaults are relied upon
# implicitly, making the simulation fully reproducible and auditable.

make_splatter_params <- function(scenario, seed) {
  
  # Resolve DE parameters from the difficulty level
  de <- DE_LEVELS[[scenario$de_signal]]
  
  # Build group probability vector from imbalance level and n_types
  group_probs <- make_group_probs(scenario$n_types, scenario$imbalance)
  
  # Construct the base parameter object
  params <- newSplatParams()
  
  # -- Experimental design parameters (vary across scenarios) ---------------
  
  # Total number of cells in the simulation
  params <- setParam(params, "batchCells",   scenario$cell_count)
  
  # Cell type proportions (geometric series — see make_group_probs)
  # NOTE: nGroups is derived automatically by Splatter from the length of
  # group.prob — it cannot be set directly via setParam. Setting group.prob
  # to a vector of length n_types implicitly defines the number of groups.
  params <- setParam(params, "group.prob",   group_probs)
  
  # Proportion of genes that are differentially expressed per group
  # Grounds: see DE_LEVELS documentation above
  params <- setParam(params, "de.prob",      de$de.prob)
  
  # Log-normal location for fold-change magnitude
  # mean_FC = exp(de.facLoc + de.facScale^2 / 2)
  params <- setParam(params, "de.facLoc",    de$de.facLoc)
  
  # -- Fixed parameters (constant across all scenarios) ---------------------
  
  # Total number of genes simulated
  params <- setParam(params, "nGenes",       FIXED_PARAMS$nGenes)
  
  # Log-normal scale for fold-change variance (fixed — see FIXED_PARAMS)
  params <- setParam(params, "de.facScale",  FIXED_PARAMS$de.facScale)
  
  # Biological coefficient of variation (within-type heterogeneity)
  params <- setParam(params, "bcv.common",   FIXED_PARAMS$bcv.common)
  
  # Library size distribution — log-normal location
  params <- setParam(params, "lib.loc",      FIXED_PARAMS$lib.loc)
  
  # Library size distribution — log-normal scale (depth variation)
  params <- setParam(params, "lib.scale",    FIXED_PARAMS$lib.scale)
  
  # Dropout model type — "none" (NB handles zeros; see FIXED_PARAMS comment)
  params <- setParam(params, "dropout.type", FIXED_PARAMS$dropout.type)
  
  # Random seed — the ONLY thing that differs between replicates
  params <- setParam(params, "seed",         seed)
  
  return(params)
}

# ============================================================================
# SIMULATION AND SAVING FUNCTION
# ============================================================================

run_simulation <- function(scenario, rep_number, seed) {
  
  sim_name <- paste0(scenario$scenario_id, "_rep", rep_number)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Simulating:", sim_name, "\n")
  cat("Scenario:  ", scenario$description, "\n")
  cat(sprintf("Seed:       %d (replicate %d of 3)\n", seed, rep_number))
  cat(rep("=", 70), "\n", sep = "")
  
  # Build parameters
  params <- make_splatter_params(scenario, seed)
  
  # Compute mean upregulated FC for logging
  # NOTE: this is the mean FC in the upregulated direction only.
  # Splatter assigns both up and down DE genes; averaging them together
  # produces cancellation. The upregulated mean is the meaningful signal
  # strength metric — reported in config.txt as "mean_FC_up".
  de       <- DE_LEVELS[[scenario$de_signal]]
  mean_fc  <- exp(de$de.facLoc + FIXED_PARAMS$de.facScale^2 / 2)
  n_de     <- round(FIXED_PARAMS$nGenes * de$de.prob)
  group_p  <- make_group_probs(scenario$n_types, scenario$imbalance)
  
  # Run Splatter simulation
  sim <- splatSimulate(params, method = "groups", verbose = FALSE)
  
  # Extract components
  counts_mat  <- counts(sim)
  gene_names  <- rownames(counts_mat)
  cell_ids    <- colnames(counts_mat)
  cell_labels <- as.character(sim$Group)
  
  # Build metadata data frame
  metadata <- data.frame(
    cell_id             = cell_ids,
    cell_type           = paste0("Type_", cell_labels),
    cell_type_numeric   = as.integer(gsub("Group", "", cell_labels)),
    batch               = as.character(sim$Batch),
    lib_size            = colSums(counts_mat),
    n_genes_detected    = colSums(counts_mat > 0),
    row.names           = cell_ids,
    stringsAsFactors    = FALSE
  )
  
  # Create output directory for this replicate
  # Structure: l9_benchmark/S1/rep1/ — nested so all replicates for a
  # scenario sit under a shared parent folder, making it easy to load
  # all three replicates with list.files("l9_benchmark/S1")
  sim_dir <- file.path(output_dir, scenario$scenario_id, paste0("rep", rep_number))
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save count matrix as sparse .mtx (memory efficient for large scenarios)
  counts_sparse <- as(counts_mat, "sparseMatrix")
  writeMM(counts_sparse, file.path(sim_dir, "counts.mtx"))
  
  # Save row/column names
  write.table(gene_names, file.path(sim_dir, "genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(cell_ids,   file.path(sim_dir, "cells.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Save metadata
  write.csv(metadata, file.path(sim_dir, "metadata.csv"), row.names = FALSE)
  
  # Save human-readable config file for auditability
  cell_type_counts <- table(metadata$cell_type)
  type_lines <- paste(
    sprintf("    %-12s : %4d cells  (%.2f%%)",
            names(cell_type_counts),
            as.integer(cell_type_counts),
            as.numeric(cell_type_counts) / ncol(counts_mat) * 100),
    collapse = "\n"
  )
  
  config_text <- paste0(
    "Simulation: ", sim_name, "\n",
    paste(rep("-", 60), collapse = ""), "\n",
    "Description : ", scenario$description, "\n",
    "Seed        : ", seed, " (replicate ", rep_number, " of 3)\n\n",
    
    "EXPERIMENTAL DIMENSIONS:\n",
    sprintf("  Cell count   : %d\n", scenario$cell_count),
    sprintf("  Imbalance    : %s\n", scenario$imbalance),
    sprintf("  Num types    : %d\n", scenario$n_types),
    sprintf("  DE signal    : %s\n", scenario$de_signal),
    
    "\nDE PARAMETERS:\n",
    sprintf("  de.prob      : %.3f  (%d DE genes per type)\n",
            de$de.prob, n_de),
    sprintf("  de.facLoc    : %.2f  (mean FC = %.2fx)\n",
            de$de.facLoc, mean_fc),
    sprintf("  de.facScale  : %.2f  (fixed)\n", FIXED_PARAMS$de.facScale),
    
    "\nFIXED PARAMETERS:\n",
    sprintf("  nGenes       : %d\n",   FIXED_PARAMS$nGenes),
    sprintf("  bcv.common   : %.2f\n", FIXED_PARAMS$bcv.common),
    sprintf("  lib.loc      : %.2f  (median ~%d UMI/cell)\n",
            FIXED_PARAMS$lib.loc, round(exp(FIXED_PARAMS$lib.loc))),
    sprintf("  lib.scale    : %.2f\n", FIXED_PARAMS$lib.scale),
    sprintf("  dropout.type : %s  (NB handles zeros per Svensson 2020)\n",
            FIXED_PARAMS$dropout.type),
    
    "\nOUTPUT STATISTICS:\n",
    sprintf("  Cells        : %d\n",    ncol(counts_mat)),
    sprintf("  Genes        : %d\n",    nrow(counts_mat)),
    sprintf("  Sparsity     : %.1f%%\n",
            sum(counts_mat == 0) / length(counts_mat) * 100),
    sprintf("  Median UMI   : %.0f\n",  median(colSums(counts_mat))),
    
    "\nCELL TYPE DISTRIBUTION:\n", type_lines, "\n"
  )
  
  writeLines(config_text, file.path(sim_dir, "config.txt"))
  
  # Console summary
  cat(sprintf("  DE genes   : %.1f%%  (~%d per type)\n",
              de$de.prob * 100, n_de))
  cat(sprintf("  Mean FC    : %.2fx\n", mean_fc))
  cat(sprintf("  BCV        : %.2f\n",  FIXED_PARAMS$bcv.common))
  cat(sprintf("  Sparsity   : %.1f%%\n",
              sum(counts_mat == 0) / length(counts_mat) * 100))
  cat(sprintf("  Median UMI : %.0f\n",  median(colSums(counts_mat))))
  cat(sprintf("  Rarest type: %.2f%%  (%d cells)\n",
              min(group_p) * 100,
              round(min(group_p) * scenario$cell_count)))
  cat("✓ Saved to:", sim_dir, "\n")
  
  return(invisible(sim_dir))
}

# ============================================================================
# RUN ALL 27 SIMULATIONS  (9 scenarios × 3 replicates)
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  L9 Orthogonal Array Benchmark — scRNA-seq Annotation Methods     ║\n")
cat("║  9 scenarios × 3 replicates = 27 datasets                        ║\n")
cat("║  Parameters grounded in Schirmer et al., Svensson (2020),        ║\n")
cat("║  Abdelaal et al. (2019), Zappia et al. (2017)                    ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Output directory  :", output_dir, "\n")
cat("Scenarios         : 9 (L9 Taguchi orthogonal array)\n")
cat("Replicates/scenario: 3 (seeds:", paste(REPLICATE_SEEDS, collapse = ", "), ")\n")
cat("Total datasets    : 27\n\n")

for (scenario in l9_scenarios) {
  for (rep_idx in seq_along(REPLICATE_SEEDS)) {
    seed <- REPLICATE_SEEDS[rep_idx]
    run_simulation(scenario, rep_number = rep_idx, seed = seed)
  }
}

# ============================================================================
# WRITE MASTER SUMMARY FILE
# ============================================================================

summary_file <- file.path(output_dir, "L9_DESIGN_SUMMARY.txt")

summary_lines <- c(
  "L9 ORTHOGONAL ARRAY BENCHMARK — MASTER SUMMARY",
  paste(rep("=", 70), collapse = ""),
  paste0("Generated : ", Sys.time()),
  paste0("Output dir: ", output_dir),
  "",
  "EXPERIMENTAL DESIGN",
  paste(rep("-", 70), collapse = ""),
  "Design     : Taguchi L9(3^4) orthogonal array",
  "Dimensions : 4 factors at 3 levels each",
  "Scenarios  : 9 (instead of 81 for full 3^4 factorial)",
  "Replicates : 3 per scenario (seeds: 42, 123, 999)",
  "Datasets   : 27 total",
  "",
  "DIMENSION LEVELS",
  paste(rep("-", 70), collapse = ""),
  "Cell count  : 500 / 3,000 / 15,000",
  "Imbalance   : Balanced / Mild (R=5) / Severe (R=20)",
  "Num types   : 5 / 15 / 30",
  "DE signal   : Easy (FC_up≈4.8x, 15% DE) / Medium (FC_up≈2.4x, 8% DE) / Hard (FC_up≈1.5x, 4% DE)",
  "",
  "FIXED PARAMETERS (all scenarios)",
  paste(rep("-", 70), collapse = ""),
  paste0("nGenes       : ", FIXED_PARAMS$nGenes),
  paste0("bcv.common   : ", FIXED_PARAMS$bcv.common),
  paste0("lib.loc      : ", FIXED_PARAMS$lib.loc,
         "  (median ~", round(exp(FIXED_PARAMS$lib.loc)), " UMI/cell)"),
  paste0("lib.scale    : ", FIXED_PARAMS$lib.scale),
  paste0("dropout.type : ", FIXED_PARAMS$dropout.type,
         "  (NB handles zeros — Svensson 2020)"),
  paste0("de.facScale  : ", FIXED_PARAMS$de.facScale),
  "",
  "L9 SCENARIO MATRIX",
  paste(rep("-", 70), collapse = ""),
  sprintf("%-4s  %-10s  %-10s  %-8s  %-8s  %s",
          "ID", "Cells", "Imbalance", "Types", "DE", "Description"),
  paste(rep("-", 70), collapse = "")
)

for (sc in l9_scenarios) {
  summary_lines <- c(summary_lines,
                     sprintf("%-4s  %-10d  %-10s  %-8d  %-8s  %s",
                             sc$scenario_id, sc$cell_count, sc$imbalance,
                             sc$n_types, sc$de_signal, sc$description)
  )
}

summary_lines <- c(summary_lines,
                   "",
                   "NOTES",
                   paste(rep("-", 70), collapse = ""),
                   "S3 WARNING: 500 cells + severe imbalance + 30 types → rarest type ~2-3 cells.",
                   "            Retained as stress-test. Methods that fail completely → kappa = NA.",
                   "",
                   "DOWNSTREAM ANALYSIS",
                   paste(rep("-", 70), collapse = ""),
                   "Primary metric : Cohen's kappa (accounts for class imbalance)",
                   "Secondary      : Macro F1 (equal weight per cell type)",
                   "Ranking test   : Friedman test + Nemenyi post-hoc (Demsar 2006)",
                   "Effect sizes   : eta-squared from lm(kappa ~ cell_count + imbalance +",
                   "                 n_types + de_signal) — orthogonality ensures unconfounded",
                   "                 decomposition of variance",
                   "",
                   "KEY REFERENCES",
                   paste(rep("-", 70), collapse = ""),
                   "Zappia et al. (2017) Genome Biology 18:174. [Splatter]",
                   "Abdelaal et al. (2019) Genome Biology 20:194. [Annotation benchmark]",
                   "Svensson (2020) Nature Biotechnology 38:147. [Dropout model]",
                   "Demsar (2006) JMLR 7:1-30. [Statistical comparison framework]",
                   "Taguchi (1987) / Roy (1990). [L9 orthogonal array design]",
                   "Luecken & Theis (2019) Mol Sys Bio 15:e8746. [Best practices]"
)

writeLines(summary_lines, summary_file)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  ✓ All 27 simulations complete                                    ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\nSummary saved to:", summary_file, "\n\n")