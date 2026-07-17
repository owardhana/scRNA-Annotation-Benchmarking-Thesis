# ============================================================================
# L9 Benchmark — Sanity Check Script
# ============================================================================
#
# PURPOSE:
#   Validates that the L9 benchmark generator produces simulated datasets
#   with the intended biological and technical properties BEFORE committing
#   to all 27 full simulation runs. Run this script first; only proceed to
#   l9_benchmark_generator.R if all checks pass.
#
# WHAT THIS CHECKS:
#   1. Sparsity          — should be 70-85% zeros (10X Chromium range)
#   2. Library size      — median UMI/cell should be 3,000-7,000
#   3. Cell proportions  — should match group.prob within sampling noise
#   4. FC distribution   — DE genes should show expected fold-change range
#                          per difficulty level (easy/medium/hard)
#   5. Visual separation — PCA plots: easy should separate cleanly, hard
#                          should show meaningful cluster overlap
#   6. S3 edge case      — flags if rarest type has < 5 cells (too few
#                          to train any classifier reliably)
#
# STRATEGY:
#   Runs 3 representative scenarios covering the full difficulty range,
#   each at replicate 1 only (seed = 42). This is sufficient to validate
#   parameters without the time cost of all 27 datasets.
#
#   S1 : Easy   — 500 cells,  balanced,  5 types  (expect clean PCA)
#   S8 : Hard   — 15000 cells, mild,     5 types  (expect overlapping PCA)
#         ^ same n_types as S1 so PCA comparison is apples-to-apples
#   S3 : Stress — 500 cells,  severe,   30 types  (expect edge case flag)
#
#   If easy (S1) and hard (S8) look similar in PCA, your de.facLoc values
#   are not creating a real difficulty gradient — stop and recalibrate.
#   If hard (S8) shows complete separation despite low FC, your bcv.common
#   is too low and clusters are unrealistically tight — increase it.
#
# EXPECTED PASS/FAIL THRESHOLDS:
#   Sparsity     : 70% - 88%   (below = too dense; above = too sparse)
#   Median UMI   : 2,000 - 8,000
#   Min type n   : >= 5 cells  (below = S3 edge case, flag WARNING not FAIL)
#   FC easy      : mean 1.8x - 2.6x  (exp(0.7 ± tolerance))
#   FC medium    : mean 1.3x - 1.8x
#   FC hard      : mean 1.1x - 1.4x
#
# OUTPUT:
#   sanity_check_plots/   — PCA and QC plots (one PDF per scenario)
#   sanity_check_report.txt — pass/fail summary
#
# NOTE:
#   This script sources the generator to reuse its parameter definitions
#   (FIXED_PARAMS, DE_LEVELS, l9_scenarios, make_splatter_params, etc.)
#   without running the full simulation loop. The source() call below
#   assumes both scripts are in the same working directory.
#
# ============================================================================

library(splatter)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Path to the main generator script — adjust if in a different directory
GENERATOR_SCRIPT <- "l9 Generator.R"

# Output directory for plots and report
CHECK_DIR <- "sanity_check_plots"
dir.create(CHECK_DIR, showWarnings = FALSE, recursive = TRUE)

# Scenarios to validate — all 9, so the gradient is checked end-to-end.
# S3 will flag a WARNING (not FAIL) for the rarest type having < 5 cells —
# this is expected by design and documented in the generator.
CHECK_SCENARIOS <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")

# Seed to use for all sanity check runs (replicate 1 only)
CHECK_SEED <- 42

# Pass/fail thresholds — updated after v3 calibration run.
# Sparsity target 72-88%: with dropout.type = "none" and lib.loc = 8.5,
# the NB distribution naturally produces ~78-82% zeros — consistent with
# the ~80% figure commonly cited for 10X Chromium data. lib.loc = 8.0
# produced 90% sparsity (too high); lib.loc = 8.5 corrects this.
# FC thresholds match recalibrated DE_LEVELS (v3):
#   easy   de.facLoc = 1.50 → mean_FC_up ≈ 4.8x
#   medium de.facLoc = 0.80 → mean_FC_up ≈ 2.4x
#   hard   de.facLoc = 0.30 → mean_FC_up ≈ 1.5x
# Reference: Luecken & Theis (2019) best practices; 10X technical notes
THRESHOLDS <- list(
  sparsity_min      = 0.72,   # below this → simulation too dense
  sparsity_max      = 0.88,   # above this → simulation too sparse
  median_umi_min    = 2000,   # below this → library sizes too low
  median_umi_max    = 10000,  # above this → library sizes too high
  min_type_cells    = 5,      # below this → edge case, flag as WARNING
  fc_easy_min       = 3.5,    # easy mean_FC_up lower bound
  fc_easy_max       = 6.5,    # easy mean_FC_up upper bound
  fc_medium_min     = 1.8,    # medium mean_FC_up lower bound
  fc_medium_max     = 3.2,    # medium mean_FC_up upper bound
  fc_hard_min       = 1.2,    # hard mean_FC_up lower bound
  fc_hard_max       = 2.0     # hard mean_FC_up upper bound
)

# ============================================================================
# LOAD GENERATOR DEFINITIONS WITHOUT RUNNING THE LOOP
# ============================================================================
#
# We parse the generator script and evaluate everything up to (but not
# including) the outer simulation loop. This gives us access to
# FIXED_PARAMS, DE_LEVELS, l9_scenarios, make_group_probs, and
# make_splatter_params without generating 27 datasets.

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  L9 Benchmark — Sanity Check                                      ║\n")
cat("║  Validating parameters before full 27-dataset run                 ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Loading generator definitions from:", GENERATOR_SCRIPT, "\n")

# Read lines of the generator, stop before the outer loop
gen_lines  <- readLines(GENERATOR_SCRIPT)
loop_start <- grep("^# RUN ALL 27 SIMULATIONS", gen_lines)[1]

if (is.na(loop_start)) {
  stop("Could not find loop header in generator script. ",
       "Check that GENERATOR_SCRIPT path is correct.")
}

# Evaluate only the definitions section
definitions_text <- paste(gen_lines[1:(loop_start - 1)], collapse = "\n")
eval(parse(text = definitions_text))
cat("✓ Definitions loaded successfully\n\n")

# ============================================================================
# HELPER: PASS / FAIL / WARNING PRINTER
# ============================================================================

check_result <- function(label, value, pass, warn = FALSE, fmt = "%.3f") {
  status <- if (pass)  "  [PASS]   "  else
    if (warn)  "  [WARN]   "  else
      "  [FAIL]   "
  cat(sprintf("%s %-40s %s\n", status, label,
              ifelse(is.numeric(value), sprintf(fmt, value), as.character(value))))
  return(pass || warn)
}

# ============================================================================
# HELPER: SINGLE-SCENARIO VALIDATION
# ============================================================================

validate_scenario <- function(scenario_id, seed = CHECK_SEED) {
  
  scenario <- l9_scenarios[[scenario_id]]
  de       <- DE_LEVELS[[scenario$de_signal]]
  
  cat("\n", rep("-", 70), "\n", sep = "")
  cat(sprintf("Scenario %s  |  %s\n", scenario_id, scenario$description))
  cat(rep("-", 70), "\n", sep = "")
  
  # -- Simulate ------------------------------------------------------------
  cat("Running splatSimulate... ")
  params <- make_splatter_params(scenario, seed)
  sim    <- splatSimulate(params, method = "groups", verbose = FALSE)
  cat("done\n\n")
  
  counts_mat  <- counts(sim)
  cell_labels <- as.character(sim$Group)
  n_cells     <- ncol(counts_mat)
  n_genes     <- nrow(counts_mat)
  group_p     <- make_group_probs(scenario$n_types, scenario$imbalance)
  
  results <- list(scenario_id = scenario_id, checks = list())
  
  # -----------------------------------------------------------------------
  # CHECK 1: Sparsity
  # -----------------------------------------------------------------------
  
  cat("CHECK 1 — Sparsity\n")
  sparsity <- sum(counts_mat == 0) / length(counts_mat)
  pass     <- sparsity >= THRESHOLDS$sparsity_min &&
    sparsity <= THRESHOLDS$sparsity_max
  
  check_result(
    label = sprintf("Sparsity (expect %.0f-%.0f%%)",
                    THRESHOLDS$sparsity_min * 100,
                    THRESHOLDS$sparsity_max * 100),
    value = sparsity * 100,
    pass  = pass,
    fmt   = "%.1f%%"
  )
  
  if (!pass && sparsity < THRESHOLDS$sparsity_min) {
    cat("         → Too dense. Lower lib.loc in FIXED_PARAMS to reduce UMI per cell\n")
    cat("           (fewer counts per cell → more NB zeros → higher sparsity)\n")
  }
  if (!pass && sparsity > THRESHOLDS$sparsity_max) {
    cat("         → Too sparse. Raise lib.loc in FIXED_PARAMS to increase UMI per cell\n")
    cat("           (more counts per cell → fewer NB zeros → lower sparsity)\n")
    cat("           Note: dropout.type = 'none' is correct — do NOT add dropout back\n")
  }
  
  results$checks$sparsity <- list(value = sparsity, pass = pass)
  
  # -----------------------------------------------------------------------
  # CHECK 2: Library size (UMI per cell)
  # -----------------------------------------------------------------------
  
  cat("\nCHECK 2 — Library Size\n")
  lib_sizes  <- colSums(counts_mat)
  median_umi <- median(lib_sizes)
  pass_umi   <- median_umi >= THRESHOLDS$median_umi_min &&
    median_umi <= THRESHOLDS$median_umi_max
  
  check_result(
    label = sprintf("Median UMI (expect %d-%d)",
                    THRESHOLDS$median_umi_min, THRESHOLDS$median_umi_max),
    value = median_umi,
    pass  = pass_umi,
    fmt   = "%.0f"
  )
  check_result(
    label = "Min UMI",
    value = min(lib_sizes),
    pass  = TRUE,   # informational only
    fmt   = "%.0f"
  )
  check_result(
    label = "Max UMI",
    value = max(lib_sizes),
    pass  = TRUE,
    fmt   = "%.0f"
  )
  
  results$checks$library_size <- list(value = median_umi, pass = pass_umi)
  
  # -----------------------------------------------------------------------
  # CHECK 3: Cell type proportions
  # -----------------------------------------------------------------------
  
  cat("\nCHECK 3 — Cell Type Proportions\n")
  obs_counts <- table(cell_labels)
  
  # CRITICAL: sort obs_counts numerically not alphabetically.
  # table() sorts group names lexicographically, which for 10+ groups
  # gives Group1, Group10, Group11, ..., Group2, Group20, ...
  # group_p is in numeric order (Group1, Group2, ..., Group30).
  # Without this sort, observed proportions for Group10 are compared
  # against Group2's expected probability, producing large false deviations
  # on any scenario with >= 10 types and unequal proportions.
  group_nums <- as.integer(gsub("Group", "", names(obs_counts)))
  obs_counts <- obs_counts[order(group_nums)]
  
  obs_props  <- as.numeric(obs_counts) / n_cells
  exp_props  <- group_p
  
  cat(sprintf("  %-12s  %8s  %8s  %8s\n",
              "Type", "Expected", "Observed", "Diff"))
  cat(sprintf("  %s\n", paste(rep("-", 44), collapse = "")))
  
  max_diff <- 0
  for (i in seq_along(obs_props)) {
    diff     <- abs(obs_props[i] - exp_props[i])
    max_diff <- max(max_diff, diff)
    cat(sprintf("  %-12s  %7.2f%%  %7.2f%%  %+7.2f%%\n",
                paste0("Group", i),
                exp_props[i] * 100,
                obs_props[i] * 100,
                (obs_props[i] - exp_props[i]) * 100))
  }
  
  # Proportions should be within ~5% of target (sampling noise)
  prop_ok <- max_diff < 0.05
  check_result(
    label = "Max deviation from target proportion",
    value = max_diff,
    pass  = prop_ok,
    fmt   = "%.3f"
  )
  
  results$checks$proportions <- list(max_diff = max_diff, pass = prop_ok)
  
  # -----------------------------------------------------------------------
  # CHECK 4: Minimum cell type count (edge case for S3)
  # -----------------------------------------------------------------------
  
  cat("\nCHECK 4 — Minimum Cell Type Count\n")
  min_cells <- min(as.integer(obs_counts))
  is_warn   <- min_cells < THRESHOLDS$min_type_cells
  is_pass   <- !is_warn
  
  check_result(
    label = sprintf("Min cells in any type (threshold: %d)",
                    THRESHOLDS$min_type_cells),
    value = min_cells,
    pass  = is_pass,
    warn  = is_warn,
    fmt   = "%.0f"
  )
  
  if (is_warn) {
    cat(sprintf("         → WARNING: Rarest type has only %d cell(s).\n",
                min_cells))
    cat("           S3 is a stress-test by design — methods may return NA.\n")
    cat("           This is expected behaviour, not a parameter error.\n")
  }
  
  results$checks$min_cells <- list(value = min_cells, pass = is_pass || is_warn)
  
  # -----------------------------------------------------------------------
  # CHECK 5: Fold-change distribution
  # -----------------------------------------------------------------------
  #
  # Splatter stores the true log-fold-change for each gene in the SCE
  # object's rowData. We extract these to verify the simulated FC values
  # match our intended de.facLoc parameter.
  
  cat("\nCHECK 5 — Fold-Change Distribution\n")
  cat("  (Validating simulated DE genes match intended difficulty level)\n\n")
  
  rd           <- rowData(sim)
  de_gene_mask <- rd$DEFacGroup1 != 1  # any gene with FC ≠ 1 is DE in Group1
  
  if (sum(de_gene_mask) == 0) {
    cat("  [WARN]    No DE genes detected in rowData — check Splatter version\n")
    results$checks$fc <- list(pass = FALSE, note = "No DE genes found")
  } else {
    # Collect all non-1 DE factors across all groups
    de_cols <- grep("^DEFac", colnames(rd), value = TRUE)
    all_fcs <- unlist(lapply(de_cols, function(col) {
      vals <- rd[[col]]
      vals[vals != 1]   # exclude non-DE genes (factor = 1)
    }))
    
    # Use only upregulated genes (FC > 1) to compute mean FC.
    # Averaging all DE genes together causes cancellation because Splatter
    # assigns both upregulated (FC > 1) and downregulated (FC < 1) genes.
    # The upregulated mean is the meaningful signal strength metric and
    # matches the de.facLoc parameter interpretation.
    up_fcs   <- all_fcs[all_fcs > 1]
    mean_fc  <- if (length(up_fcs) > 0) mean(up_fcs) else NA
    q25_fc   <- if (length(up_fcs) > 0) quantile(up_fcs, 0.25) else NA
    q75_fc   <- if (length(up_fcs) > 0) quantile(up_fcs, 0.75) else NA
    n_de_obs <- length(all_fcs)
    
    # Determine expected range from the scenario's DE level
    fc_min <- THRESHOLDS[[paste0("fc_", scenario$de_signal, "_min")]]
    fc_max <- THRESHOLDS[[paste0("fc_", scenario$de_signal, "_max")]]
    pass_fc <- mean_fc >= fc_min && mean_fc <= fc_max
    
    cat(sprintf("  DE level  : %s  (target mean_FC_up: %.1fx - %.1fx)\n",
                scenario$de_signal, fc_min, fc_max))
    cat(sprintf("  n DE genes: %d  (%d upregulated)\n", n_de_obs, length(up_fcs)))
    
    check_result(
      label = sprintf("Mean FC upregulated (expect %.1fx - %.1fx)", fc_min, fc_max),
      value = mean_fc,
      pass  = pass_fc,
      fmt   = "%.2fx"
    )
    check_result(
      label = "25th percentile FC (upregulated)",
      value = q25_fc,
      pass  = TRUE,
      fmt   = "%.2fx"
    )
    check_result(
      label = "75th percentile FC (upregulated)",
      value = q75_fc,
      pass  = TRUE,
      fmt   = "%.2fx"
    )
    
    if (!pass_fc && !is.na(mean_fc) && mean_fc < fc_min) {
      cat(sprintf("         → FC too low. de.facLoc = %.2f may need increasing.\n",
                  de$de.facLoc))
    }
    if (!pass_fc && !is.na(mean_fc) && mean_fc > fc_max) {
      cat(sprintf("         → FC too high. de.facLoc = %.2f may need decreasing.\n",
                  de$de.facLoc))
    }
    
    results$checks$fc <- list(mean_fc = mean_fc, pass = pass_fc)
  }
  
  # -----------------------------------------------------------------------
  # CHECK 6: PCA visual separation
  # -----------------------------------------------------------------------
  #
  # Runs PCA on log-normalized counts, plots PC1 vs PC2 coloured by
  # cell type. For easy scenarios, types should separate visibly in
  # the top 2 PCs. For hard scenarios, substantial overlap is expected
  # and desired — if hard looks as clean as easy, the difficulty gradient
  # is not working.
  #
  # This check is VISUAL ONLY — no numeric pass/fail threshold, since
  # "separation" is context-dependent. The plot is saved to CHECK_DIR
  # and must be inspected manually.
  
  cat("\nCHECK 6 — PCA Visual Separation (inspect plot manually)\n")
  
  # Log-normalize counts
  lib_size_factors <- colSums(counts_mat) / median(colSums(counts_mat))
  norm_counts      <- log1p(t(t(counts_mat) / lib_size_factors))
  
  # Select top 2000 most variable genes for PCA efficiency
  gene_vars    <- apply(norm_counts, 1, var)
  top_genes    <- order(gene_vars, decreasing = TRUE)[1:min(2000, nrow(norm_counts))]
  pca_input    <- t(norm_counts[top_genes, ])
  
  # Run PCA
  pca_result   <- prcomp(pca_input, center = TRUE, scale. = FALSE)
  pca_df       <- data.frame(
    PC1       = pca_result$x[, 1],
    PC2       = pca_result$x[, 2],
    CellType  = factor(cell_labels),
    LibSize   = lib_sizes
  )
  
  # Variance explained
  var_exp <- summary(pca_result)$importance[2, 1:min(10, ncol(pca_result$x))]
  
  cat(sprintf("  PC1 variance explained: %.1f%%\n", var_exp[1] * 100))
  cat(sprintf("  PC2 variance explained: %.1f%%\n", var_exp[2] * 100))
  cat(sprintf("  PC1+PC2 combined      : %.1f%%\n",
              sum(var_exp[1:2]) * 100))
  
  if (sum(var_exp[1:2]) < 0.10) {
    cat("  [WARN]    < 10% variance in PC1+PC2. Clusters likely diffuse —\n")
    cat("            check BCV and DE parameters if this is the easy scenario.\n")
  }
  
  # -- Plot 1: PCA coloured by cell type
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = CellType)) +
    geom_point(alpha = 0.6, size = 1.2) +
    labs(
      title    = sprintf("PCA — Scenario %s (%s DE)", scenario_id, scenario$de_signal),
      subtitle = sprintf("%d cells, %d types, %s imbalance | PC1+PC2: %.1f%% var",
                         n_cells, scenario$n_types, scenario$imbalance,
                         sum(var_exp[1:2]) * 100),
      x        = sprintf("PC1 (%.1f%%)", var_exp[1] * 100),
      y        = sprintf("PC2 (%.1f%%)", var_exp[2] * 100),
      colour   = "Cell Type"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
  
  # -- Plot 2: Library size distribution
  lib_df <- data.frame(LibSize = lib_sizes, CellType = factor(cell_labels))
  
  p_lib <- ggplot(lib_df, aes(x = LibSize)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = median(lib_sizes),
               linetype = "dashed", colour = "red", linewidth = 0.8) +
    annotate("text",
             x     = median(lib_sizes) * 1.05,
             y     = Inf, vjust = 2, hjust = 0,
             label = sprintf("Median: %.0f", median(lib_sizes)),
             size  = 3.5, colour = "red") +
    labs(
      title    = sprintf("Library Size Distribution — %s", scenario_id),
      subtitle = sprintf("Sparsity: %.1f%%  |  Median UMI: %.0f",
                         sparsity * 100, median_umi),
      x        = "UMI per Cell",
      y        = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
  
  # -- Plot 3: Cell type proportions (observed vs expected)
  prop_df <- data.frame(
    Type     = rep(paste0("Group", seq_along(obs_props)), 2),
    Prop     = c(obs_props, exp_props) * 100,
    Source   = rep(c("Observed", "Expected"), each = length(obs_props))
  )
  
  p_prop <- ggplot(prop_df, aes(x = Type, y = Prop, fill = Source)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = c("Observed" = "steelblue",
                                 "Expected" = "coral")) +
    labs(
      title    = sprintf("Cell Type Proportions — %s", scenario_id),
      subtitle = sprintf("Max deviation: %.2f%% | Min type: %d cells",
                         max_diff * 100, min_cells),
      x        = "Cell Type",
      y        = "Proportion (%)",
      fill     = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          plot.title   = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
  
  # -- Combine and save
  combined_plot <- (p_pca | p_lib) / p_prop +
    plot_annotation(
      title   = sprintf("Sanity Check — Scenario %s", scenario_id),
      caption = "Red dashed line = median. Inspect PCA for expected separation level.",
      theme   = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  plot_file <- file.path(CHECK_DIR, sprintf("sanity_%s.pdf", scenario_id))
  ggsave(plot_file, combined_plot, width = 14, height = 10)
  cat(sprintf("  Plot saved to: %s\n", plot_file))
  
  # -- Visual guidance for interpreter
  cat("\n  PCA INTERPRETATION GUIDE:\n")
  if (scenario$de_signal == "easy") {
    cat("  This is an EASY scenario — cell types should form clearly\n")
    cat("  separated clusters in PC1/PC2. If they overlap substantially,\n")
    cat("  de.facLoc = 0.70 is not producing enough signal — recalibrate.\n")
  } else if (scenario$de_signal == "hard") {
    cat("  This is a HARD scenario — some cluster overlap is EXPECTED and\n")
    cat("  DESIRED. If clusters look as clean as the easy scenario (S1),\n")
    cat("  de.facLoc = 0.12 is not producing sufficient difficulty.\n")
    cat("  Try reducing de.facLoc further or increasing bcv.common.\n")
  }
  
  results$checks$pca <- list(
    pc1_var = var_exp[1],
    pc2_var = var_exp[2],
    pass    = TRUE  # visual check — always marked pass, review plot manually
  )
  
  return(results)
}

# ============================================================================
# RUN CHECKS ON THE THREE REPRESENTATIVE SCENARIOS
# ============================================================================

cat("Checking scenarios:", paste(CHECK_SCENARIOS, collapse = ", "), "\n")
cat("Seed used:", CHECK_SEED, "(replicate 1 only)\n")

all_results <- list()

for (sid in CHECK_SCENARIOS) {
  all_results[[sid]] <- validate_scenario(sid, seed = CHECK_SEED)
}

# ============================================================================
# CROSS-SCENARIO COMPARISON: FC GRADIENT CHECK
# ============================================================================
#
# The most important validation: does the difficulty gradient actually
# produce meaningfully different levels of cluster separability?
# We compare PC1 variance explained across S1 (easy) and S8 (hard).
# Easy should have substantially more variance in PC1 than hard.

cat("\n", rep("=", 70), "\n", sep = "")
cat("CROSS-SCENARIO COMPARISON — Difficulty Gradient\n")
cat(rep("=", 70), "\n", sep = "")

# Compare PC1 variance across the three DE signal levels.
# S1 = easy, S8 = hard — both have n_types=5 so PCA is apples-to-apples.
# S5 = easy with 30 types (different structure, for reference only).
# The key question: does easy have substantially more PC1 variance than hard?

s1_pc1 <- all_results[["S1"]]$checks$pca$pc1_var
s8_pc1 <- all_results[["S8"]]$checks$pca$pc1_var

# Also extract medium-DE scenario with 5 types for a full gradient check.
# S6 has 5 types + medium DE — most comparable to S1 (easy) and S8 (hard).
s6_pc1 <- if ("S6" %in% names(all_results)) {
  all_results[["S6"]]$checks$pca$pc1_var
} else NA

cat(sprintf("\n  PC1 variance explained by DE level (5-type scenarios):\n"))
cat(sprintf("    S1 (easy,   5 types, balanced) : %.1f%%\n", s1_pc1 * 100))
if (!is.na(s6_pc1)) {
  cat(sprintf("    S6 (medium, 5 types, severe)   : %.1f%%\n", s6_pc1 * 100))
}
cat(sprintf("    S8 (hard,   5 types, mild)     : %.1f%%\n", s8_pc1 * 100))
cat(sprintf("    Ratio easy/hard: %.1fx\n", s1_pc1 / s8_pc1))

if (s1_pc1 > s8_pc1 * 1.5) {
  cat("\n  [PASS]    Easy scenario has substantially more PC1 variance than hard.\n")
  cat("            Difficulty gradient is working as intended.\n")
} else if (s1_pc1 > s8_pc1 * 1.1) {
  cat("\n  [WARN]    Easy scenario has only marginally more PC1 variance than hard.\n")
  cat("            Gradient may be too weak. Consider widening de.facLoc gap.\n")
} else {
  cat("\n  [FAIL]    Easy and hard scenarios have similar PC1 variance.\n")
  cat("            The difficulty gradient is NOT working. Do not proceed\n")
  cat("            with the full 27-dataset run until parameters are fixed.\n")
}

# ============================================================================
# SUMMARY REPORT
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("OVERALL SANITY CHECK SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

report_lines <- c(
  "L9 BENCHMARK — SANITY CHECK REPORT",
  paste(rep("=", 70), collapse = ""),
  paste0("Generated : ", Sys.time()),
  paste0("Generator : ", GENERATOR_SCRIPT),
  paste0("Scenarios : ", paste(CHECK_SCENARIOS, collapse = ", ")),
  paste0("Seed      : ", CHECK_SEED, " (replicate 1 only)"),
  ""
)

all_passed <- TRUE

for (sid in CHECK_SCENARIOS) {
  r      <- all_results[[sid]]
  checks <- r$checks
  
  scenario_pass <- all(sapply(checks, function(x) isTRUE(x$pass)))
  all_passed    <- all_passed && scenario_pass
  
  status_str <- if (scenario_pass) "[PASS]" else "[FAIL/WARN — review above]"
  
  cat(sprintf("  %s  %s\n", status_str, sid))
  report_lines <- c(report_lines,
                    sprintf("%s  %s", status_str, sid),
                    sprintf("    Sparsity     : %.1f%%  (threshold: %.0f-%.0f%%)",
                            checks$sparsity$value * 100,
                            THRESHOLDS$sparsity_min * 100, THRESHOLDS$sparsity_max * 100),
                    sprintf("    Median UMI   : %.0f  (threshold: %d-%d)",
                            checks$library_size$value,
                            THRESHOLDS$median_umi_min, THRESHOLDS$median_umi_max),
                    sprintf("    Min type n   : %d  (threshold: >=%d)",
                            checks$min_cells$value, THRESHOLDS$min_type_cells),
                    ifelse(!is.null(checks$fc$mean_fc),
                           sprintf("    Mean FC_up   : %.2fx  (upregulated direction only)", checks$fc$mean_fc),
                           "    Mean FC_up   : not computed"),
                    ""
  )
}

report_lines <- c(report_lines,
                  "",
                  "DIFFICULTY GRADIENT CHECK (5-type scenarios: S1=easy, S6=medium, S8=hard)",
                  paste(rep("-", 70), collapse = ""),
                  sprintf("S1 (easy)   PC1 var : %.1f%%", s1_pc1 * 100),
                  sprintf("S6 (medium) PC1 var : %.1f%%", ifelse(is.na(s6_pc1), NA, s6_pc1 * 100)),
                  sprintf("S8 (hard)   PC1 var : %.1f%%", s8_pc1 * 100),
                  sprintf("Ratio easy/hard     : %.1fx",  s1_pc1 / s8_pc1),
                  "",
                  "PLOTS (inspect manually — all 9 scenarios)",
                  paste(rep("-", 70), collapse = "")
)

for (sid in CHECK_SCENARIOS) {
  report_lines <- c(report_lines,
                    sprintf("  %s/sanity_%s.pdf", CHECK_DIR, sid)
  )
}

report_lines <- c(report_lines,
                  "",
                  "NEXT STEPS",
                  paste(rep("-", 70), collapse = ""),
                  "If all checks pass AND PCA plots look correct:",
                  "  → Run l9_benchmark_generator.R to generate all 27 datasets",
                  "",
                  "If sparsity/UMI checks fail:",
                  "  → Adjust dropout.mid or lib.loc in FIXED_PARAMS",
                  "",
                  "If FC checks fail:",
                  "  → Adjust de.facLoc in DE_LEVELS",
                  "",
                  "If difficulty gradient is flat (easy ≈ hard in PCA):",
                  "  → Widen the de.facLoc gap between easy and hard levels",
                  "  → Or increase bcv.common to add realistic within-type noise"
)

report_file <- "sanity_check_report.txt"
writeLines(report_lines, report_file)

cat("\n")
if (all_passed) {
  cat("╔════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✓ All numeric checks passed                                      ║\n")
  cat("║  Review PCA plots manually before running full generator          ║\n")
  cat("╚════════════════════════════════════════════════════════════════════╝\n")
} else {
  cat("╔════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✗ One or more checks failed — do NOT run full generator yet      ║\n")
  cat("║  Review output above and adjust parameters in generator script    ║\n")
  cat("╚════════════════════════════════════════════════════════════════════╝\n")
}

cat(sprintf("\nReport saved to : %s\n", report_file))
cat(sprintf("Plots saved to  : %s/\n\n", CHECK_DIR))