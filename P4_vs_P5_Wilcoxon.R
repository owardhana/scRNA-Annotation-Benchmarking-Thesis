# ============================================================
# P4 (Synthetic) vs P5 V2 (Real Life) — Wilcoxon Signed-Rank Test
# Paired comparison of benchmark metrics per tool category
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)

# ============================================================
# 1. CONSTANTS
# ============================================================

CATEGORY_COLOURS <- c(
  "Marker-Based"      = "#E07B39",
  "Correlation-Based"  = "#4A90C4",
  "Classic ML"         = "#5BAD6F",
  "Deep Learning"      = "#9B6BB5",
  "Semi-Supervised"    = "#1A9C8E"
)

EXCLUDE_TOOLS <- c(
  'scInfeR', 'scPred_adaboost', "scAnnotatR", "CIPR",
  "scmap_cluster", "clustifyr_hyper", "clustifyr_jaccard",
  "SCSA", "scCATCH", "scType", "scGAD", "scDeepSort", "CALLR"
)

PLOT_DIR <- "./Plots"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

# ============================================================
# 2. TOOL CATEGORY DEFINITIONS
# ============================================================

marker_tools <- c("clustifyr_hyper","clustifyr_jaccard","scCATCH","SCINA",
                  "ScInfeR","SCSA","scSorter","scType")

correlation_tools <- c("CIPR","scibetR","scmap_cell","scmap_cluster",
                       "Seurat_Transfer_CCA","Seurat_Transfer_PCA",
                       "Seurat_Transfer_RPCA","SingleR")

classicML_tools <- c("CaSTLe","CellTypist","CHETAH","scAnnotate","scAnnotatR",
                     "scClassify","scID","scPred","scPred_adaboost","scPred_avNNet",
                     "scPred_bayesglm","scPred_earth","scPred_glm","scPred_glmboost",
                     "scPred_knn","scPred_lda","scPred_mlp","scPred_multinom",
                     "scPred_nb","scPred_nnet","scPred_regLogistic","scPred_rf",
                     "scPred_svmLinear","scPred_svmPoly","scPred_xgbTree",
                     "singleCellNet")

DL_tools <- c("ACTINN","Cell_BLAST","CIForm","NeuCA","scBalance","scDeepSort",
              "scHash","scLearn","scMMT","TOSICA")

semi_supervised_tools <- c("CALLR","CAMLU","MARS","SCANVI","scArches","scGAD",
                           "scnym","scSemiCluster","scSemiGAN")

.build_category_map <- function() {
  c(setNames(rep("Marker-Based",      length(marker_tools)),          marker_tools),
    setNames(rep("Correlation-Based",  length(correlation_tools)),     correlation_tools),
    setNames(rep("Classic ML",         length(classicML_tools)),       classicML_tools),
    setNames(rep("Deep Learning",      length(DL_tools)),              DL_tools),
    setNames(rep("Semi-Supervised",    length(semi_supervised_tools)), semi_supervised_tools))
}

.drop_index_col <- function(df) {
  if (ncol(df) > 0 && colnames(df)[1] %in% c("", "X", "Unnamed: 0"))
    df[, -1, drop = FALSE] else df
}

# ============================================================
# 3. DATA LOADER (supports both nested-replicate & flat layouts)
# ============================================================

load_benchmark_results <- function(base_dir, subfolders = NULL,
                                   exclude_tools = NULL, pattern = "\\.csv$") {
  if (!dir.exists(base_dir)) stop("base_dir does not exist: ", base_dir)
  if (is.null(subfolders)) {
    subfolders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
    subfolders <- subfolders[nchar(subfolders) > 0]
  }
  cat_map  <- .build_category_map()
  all_dfs  <- list()

  for (scenario in subfolders) {
    sp <- file.path(base_dir, scenario)
    if (!dir.exists(sp)) next

    rep_dirs <- list.dirs(sp, full.names = FALSE, recursive = FALSE)
    rep_dirs <- rep_dirs[nchar(rep_dirs) > 0]

    if (length(rep_dirs) == 0) {
      # Flat layout (P5 V2 style — CSVs directly in scenario folder)
      csvs <- list.files(sp, pattern = pattern, full.names = TRUE, recursive = FALSE)
      if (length(csvs) == 0) next
      dfs <- lapply(csvs, function(f) tryCatch({
        df <- .drop_index_col(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
        df$scenario  <- scenario
        df$replicate <- "rep1"
        df
      }, error = function(e) NULL))
      dfs <- Filter(Negate(is.null), dfs)
      if (length(dfs) > 0) all_dfs[[scenario]] <- bind_rows(dfs)
    } else {
      # Nested replicate layout (P4 style)
      for (rep in rep_dirs) {
        rp   <- file.path(sp, rep)
        csvs <- list.files(rp, pattern = pattern, full.names = TRUE, recursive = FALSE)
        if (length(csvs) == 0) next
        dfs <- lapply(csvs, function(f) tryCatch({
          df <- .drop_index_col(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
          df$scenario  <- scenario
          df$replicate <- rep
          df
        }, error = function(e) NULL))
        dfs <- Filter(Negate(is.null), dfs)
        if (length(dfs) > 0) all_dfs[[paste0(scenario, "/", rep)]] <- bind_rows(dfs)
      }
    }
  }

  if (length(all_dfs) == 0) stop("No data loaded.")
  combined <- bind_rows(all_dfs)
  rownames(combined) <- NULL
  if (!is.null(exclude_tools))
    combined <- combined[!combined$tool %in% exclude_tools, , drop = FALSE]
  combined$tool_category <- cat_map[combined$tool]
  combined
}

# ============================================================
# 4. LOAD P4 & P5 V2
# ============================================================

P4_DIR <- "../P4 Synthetic Datasets V1"
P5_DIR <- "./"

df_p4 <- load_benchmark_results(
  base_dir      = P4_DIR,
  subfolders    = paste0("S", 1:9),
  exclude_tools = EXCLUDE_TOOLS
)

message("P4 loaded: ", nrow(df_p4), " rows, columns: ", paste(colnames(df_p4), collapse=", "))

df_p5 <- load_benchmark_results(
  base_dir      = P5_DIR,
  subfolders    = c(
    "S1-Darmanis-Brain-2015_for_use",
    "S2-Marques-Brain_for_use",
    "S3-Nowakowski-Cortex-2017_for_use",
    "S4-Grun-Pancreas-2016_for_use",
    "S5-Tabula Muris-FACS-3k_for_use",
    "S6-He-Skin-2020_for_use",
    "S7-Zhao-Immune-Fine-2020_for_use",
    "S8-Zheng-Zhengsort-5cl-2017_original",
    "S9-MacParlan-Liver-Broad_original"
  ),
  exclude_tools = EXCLUDE_TOOLS
)

message("P5 loaded: ", nrow(df_p5), " rows, columns: ", paste(colnames(df_p5), collapse=", "))

message("P4 tools: ", length(unique(df_p4$tool)),
        " | P5 tools: ", length(unique(df_p5$tool)))

# ============================================================
# 5. COMPUTE TOOL-LEVEL MEANS (one row per tool per pipeline)
# ============================================================

summarise_pipeline <- function(df, pipeline_name) {
  if (nrow(df) == 0) {
    warning("Empty dataframe for ", pipeline_name)
    return(tibble(tool=character(), tool_category=character(), kappa_mean=numeric(), f1_mean=numeric(), runtime_mean=numeric(), memory_mean=numeric(), pipeline=character()))
  }
  
  # Ensure necessary columns exist (fill with NA if missing)
  for (col in c("kappa_mean", "f1_mean", "runtime_mean", "peak_system_memory_mean_mb")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA_real_
      warning("Column ", col, " missing in ", pipeline_name, ", filling with NA")
    }
  }

  df %>%
    dplyr::filter(!is.na(tool_category)) %>%
    dplyr::group_by(tool, tool_category) %>%
    dplyr::summarise(
      kappa_mean   = mean(kappa_mean,   na.rm = TRUE),
      f1_mean      = mean(f1_mean,      na.rm = TRUE),
      runtime_mean = mean(runtime_mean, na.rm = TRUE),
      memory_mean  = mean(peak_system_memory_mean_mb, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(pipeline = pipeline_name)
}

means_p4 <- summarise_pipeline(df_p4, "Synthetic (P4)")
means_p5 <- summarise_pipeline(df_p5, "Real Life (P5)")

# Keep only tools present in BOTH pipelines (required for paired test)
shared_tools <- intersect(means_p4$tool, means_p5$tool)
message("Shared tools for pairing: ", length(shared_tools))

paired <- dplyr::inner_join(
  means_p4 %>% dplyr::filter(tool %in% shared_tools) %>%
    dplyr::rename(kappa_P4 = kappa_mean, f1_P4 = f1_mean,
                  runtime_P4 = runtime_mean, memory_P4 = memory_mean) %>%
    dplyr::select(tool, tool_category, kappa_P4, f1_P4, runtime_P4, memory_P4),
  means_p5 %>% dplyr::filter(tool %in% shared_tools) %>%
    dplyr::rename(kappa_P5 = kappa_mean, f1_P5 = f1_mean,
                  runtime_P5 = runtime_mean, memory_P5 = memory_mean) %>%
    dplyr::select(tool, kappa_P5, f1_P5, runtime_P5, memory_P5),
  by = "tool"
)

# ============================================================
# 6. WILCOXON SIGNED-RANK TEST (per category × metric)
# ============================================================

METRICS <- c("kappa", "f1", "runtime", "memory")
METRIC_LABELS <- c(kappa   = "Cohen's \u03ba",
                   f1      = "Macro F1",
                   runtime = "Runtime (s)",
                   memory  = "Peak Memory (MB)")

run_wilcoxon <- function(data, metric) {
  p4_col <- paste0(metric, "_P4")
  p5_col <- paste0(metric, "_P5")

  results <- data %>%
    dplyr::group_by(tool_category) %>%
    dplyr::group_modify(function(grp, key) {
      x <- grp[[p4_col]]
      y <- grp[[p5_col]]
      # Remove NA pairs
      valid <- complete.cases(x, y)
      x <- x[valid]; y <- y[valid]
      n <- length(x)
      if (n < 3) {
        return(tibble(
          n_tools = n, median_P4 = NA_real_, median_P5 = NA_real_,
          median_diff = NA_real_, p_value = NA_real_,
          statistic = NA_real_, direction = "insufficient data"
        ))
      }
      wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
      tibble(
        n_tools     = n,
        median_P4   = median(x, na.rm = TRUE),
        median_P5   = median(y, na.rm = TRUE),
        median_diff = median(x - y, na.rm = TRUE),
        p_value     = wt$p.value,
        statistic   = wt$statistic,
        direction   = ifelse(median(x) > median(y), "P4 > P5", "P5 > P4")
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      metric   = metric,
      sig_star = dplyr::case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  results
}

wilcoxon_results <- bind_rows(lapply(METRICS, function(m) run_wilcoxon(paired, m)))

# Also run a pooled (all categories combined) test per metric
wilcoxon_pooled <- bind_rows(lapply(METRICS, function(m) {
  p4_col <- paste0(m, "_P4")
  p5_col <- paste0(m, "_P5")
  x <- paired[[p4_col]]; y <- paired[[p5_col]]
  valid <- complete.cases(x, y); x <- x[valid]; y <- y[valid]
  if (length(x) < 3) return(NULL)
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
  tibble(
    tool_category = "ALL (pooled)", metric = m,
    n_tools = length(x),
    median_P4 = median(x), median_P5 = median(y),
    median_diff = median(x - y),
    p_value = wt$p.value, statistic = wt$statistic,
    direction = ifelse(median(x) > median(y), "P4 > P5", "P5 > P4"),
    sig_star = case_when(
      wt$p.value < 0.001 ~ "***", wt$p.value < 0.01  ~ "**",
      wt$p.value < 0.05  ~ "*",   TRUE ~ "ns"
    )
  )
}))

all_results <- bind_rows(wilcoxon_results, wilcoxon_pooled)

# ============================================================
# 7. CONSOLE OUTPUT
# ============================================================

message("\n", paste(rep("\u2550", 70), collapse = ""))
message("  WILCOXON SIGNED-RANK TEST: P4 (Synthetic) vs P5 V2 (Real Life)")
message(paste(rep("\u2550", 70), collapse = ""))

for (m in METRICS) {
  subset <- all_results %>% filter(metric == m)
  message(sprintf("\n  ── %s ──", METRIC_LABELS[m]))
  for (i in seq_len(nrow(subset))) {
    r <- subset[i, ]
    if (r$direction == "insufficient data") {
      message(sprintf("    %-20s  n=%d  (insufficient data for test)",
                      r$tool_category, r$n_tools))
    } else {
      message(sprintf(
        "    %-20s  n=%-3d  P4 med=%.3f  P5 med=%.3f  diff=%.3f  W=%.0f  p=%.4f %s  [%s]",
        r$tool_category, r$n_tools, r$median_P4, r$median_P5,
        r$median_diff, r$statistic, r$p_value, r$sig_star, r$direction
      ))
    }
  }
}

message("\n", paste(rep("\u2550", 70), collapse = ""))

# ============================================================
# 8. PLOTS
# ============================================================

theme_bench <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13, colour = "grey10"),
      plot.subtitle = element_text(colour = "grey40", size = 10, margin = margin(b = 8)),
      strip.text    = element_text(face = "bold", size = 10),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# --- 8a. Paired difference dot plot (kappa and F1) ---
paired_long <- paired %>%
  pivot_longer(
    cols      = c(kappa_P4, kappa_P5),
    names_to  = c("metric", "pipeline"),
    names_pattern = "(kappa)_(P4|P5)",
    values_to = "value"
  ) %>%
  mutate(
    pipeline = recode(pipeline, P4 = "Synthetic (P4)", P5 = "Real Life (P5)"),
    metric   = recode(metric, kappa = "Cohen's \u03ba"),
    tool_category = factor(tool_category, levels = names(CATEGORY_COLOURS))
  )

p_paired <- paired_long %>%
  ggplot(aes(x = value, y = reorder(tool, value), colour = pipeline)) +
  geom_line(aes(group = interaction(tool, metric)),
            colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.85) +
  facet_grid(tool_category ~ metric, scales = "free", space = "free_y") +
  scale_colour_manual(
    values = c("Synthetic (P4)" = "#3498DB", "Real Life (P5)" = "#E74C3C"),
    name   = "Pipeline"
  ) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Paired tool performance: Synthetic (P4) vs Real Life (P5 V2)",
    subtitle = "Each line connects the same tool across pipelines. Faceted by category and metric.",
    x = "Score", y = NULL
  ) +
  theme_bench() +
  theme(axis.text.y = element_text(size = 7),
        strip.text.y = element_text(angle = 0, hjust = 0))

# --- 8b. Category-level summary data for significance stars ---
diff_data <- paired %>%
  mutate(
    kappa_diff = kappa_P4 - kappa_P5,
    tool_category = factor(tool_category, levels = names(CATEGORY_COLOURS))
  ) %>%
  pivot_longer(
    cols      = c(kappa_diff),
    names_to  = "metric",
    values_to = "diff"
  ) %>%
  mutate(metric = recode(metric,
                         kappa_diff = "Cohen's \u03ba difference"))

max_diff_vals <- diff_data %>%
  group_by(metric) %>%
  summarise(max_val = max(diff, na.rm = TRUE), .groups = "drop")

sig_bars <- wilcoxon_results %>%
  filter(!is.na(p_value), metric %in% c("kappa")) %>%
  mutate(
    metric = recode(metric, kappa = "Cohen's \u03ba difference"),
    tool_category = factor(tool_category, levels = names(CATEGORY_COLOURS))
  ) %>%
  left_join(max_diff_vals, by = "metric")

# --- 8c. Boxplot of per-tool differences by category (with significance stars) ---
p_box <- diff_data %>%
  ggplot(aes(x = tool_category, y = diff,
             fill = tool_category, colour = tool_category)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.7) +
  geom_text(data = sig_bars, aes(x = tool_category, y = max_val + 0.05, label = sig_star),
            colour = "grey20", size = 5, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~metric) +
  scale_fill_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  labs(
    title    = "Per-tool score change: Synthetic (P4) \u2212 Real Life (P5 V2)",
    subtitle = "Each point = one tool. Values above zero = better on synthetic data. Stars = Wilcoxon signed-rank p-value.",
    x = NULL, y = "Score difference (P4 \u2212 P5)"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Render and save ---
plot(p_paired)
plot(p_box)

ggsave("P4vP5_wilcoxon_paired_dotplot.png",    p_paired,  width = 14, height = 16, dpi = 180, path = PLOT_DIR)
ggsave("P4vP5_wilcoxon_diff_boxplot.png",      p_box,     width = 12, height = 6,  dpi = 180, path = PLOT_DIR)

message("\nAll plots saved to: ", PLOT_DIR)

# ============================================================
# 9. EXPORT RESULTS TABLE
# ============================================================

write.csv(all_results, file.path(PLOT_DIR, "P4vP5_wilcoxon_results_table.csv"), row.names = FALSE)
message("Results table saved to: ", file.path(PLOT_DIR, "wilcoxon_results_table.csv"))
