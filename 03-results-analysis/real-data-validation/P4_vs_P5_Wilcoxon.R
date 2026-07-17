# ============================================================
# P1 (Synthetic) vs P2 V2 (Real Life) — Wilcoxon Signed-Rank Test
# Paired comparison of benchmark metrics per tool category
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)
library(extrafont)   # Arial font support for cairo_pdf
library(Cairo)       # cairo_pdf device

# ============================================================
# 1. CONSTANTS
# ============================================================

CATEGORY_COLOURS <- c(
  "Marker-Based"      = "#E69F00",
  "Correlation-Based"  = "#56B4E9",
  "Classic ML"         = "#009E73",
  "Deep Learning"      = "#CC79A7",
  "Semi-Supervised"    = "#0072B2"
)

EXCLUDE_TOOLS <- c(
  'ScInfeR', 'scPred_adaboost', "scAnnotatR", "CIPR",
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
      # Flat layout (P2 V2 style — CSVs directly in scenario folder)
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
      # Nested replicate layout (P1 style)
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
# 4. LOAD P1 & P2 V2
# ============================================================

P1_DIR <- "../P4 Synthetic Datasets V1"
P2_DIR <- "./"

df_P1 <- load_benchmark_results(
  base_dir      = P1_DIR,
  subfolders    = paste0("S", 1:9),
  exclude_tools = EXCLUDE_TOOLS
)

message("P1 loaded: ", nrow(df_P1), " rows, columns: ", paste(colnames(df_P1), collapse=", "))

df_P2 <- load_benchmark_results(
  base_dir      = P2_DIR,
  subfolders    = c(
    "S1-Darmanis-Brain-2015_for_use",
    "S2-Marques-Brain_for_use",
    "S3-Nowakowski-Cortex-2017_for_use",
    "S4-Grun-Pancreas-2016_original",
    "S5-Tabula Muris-FACS-3k_for_use",
    "S6-He-Skin-2020_for_use",
    "S7-Zhao-Immune-Fine-2020_for_use",
    "S8-Zheng-Zhengsort-5cl-2017_original",
    "S9-MacParland-Liver-Broad_original"
  ),
  exclude_tools = EXCLUDE_TOOLS
)

message("P2 loaded: ", nrow(df_P2), " rows, columns: ", paste(colnames(df_P2), collapse=", "))

message("P1 tools: ", length(unique(df_P1$tool)),
        " | P2 tools: ", length(unique(df_P2$tool)))

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

means_P1 <- summarise_pipeline(df_P1, "Synthetic (P1)")
means_P2 <- summarise_pipeline(df_P2, "Real Life (P2)")

# Keep only tools present in BOTH pipelines (required for paired test)
shared_tools <- intersect(means_P1$tool, means_P2$tool)
message("Shared tools for pairing: ", length(shared_tools))

paired <- dplyr::inner_join(
  means_P1 %>% dplyr::filter(tool %in% shared_tools) %>%
    dplyr::rename(kappa_P1 = kappa_mean, f1_P1 = f1_mean,
                  runtime_P1 = runtime_mean, memory_P1 = memory_mean) %>%
    dplyr::select(tool, tool_category, kappa_P1, f1_P1, runtime_P1, memory_P1),
  means_P2 %>% dplyr::filter(tool %in% shared_tools) %>%
    dplyr::rename(kappa_P2 = kappa_mean, f1_P2 = f1_mean,
                  runtime_P2 = runtime_mean, memory_P2 = memory_mean) %>%
    dplyr::select(tool, kappa_P2, f1_P2, runtime_P2, memory_P2),
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
  P1_col <- paste0(metric, "_P1")
  P2_col <- paste0(metric, "_P2")

  results <- data %>%
    dplyr::group_by(tool_category) %>%
    dplyr::group_modify(function(grp, key) {
      x <- grp[[P1_col]]
      y <- grp[[P2_col]]
      # Remove NA pairs
      valid <- complete.cases(x, y)
      x <- x[valid]; y <- y[valid]
      n <- length(x)
      if (n < 3) {
        return(tibble(
          n_tools = n, median_P1 = NA_real_, median_P2 = NA_real_,
          median_diff = NA_real_, p_value = NA_real_,
          statistic = NA_real_, direction = "insufficient data"
        ))
      }
      wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
      tibble(
        n_tools     = n,
        median_P1   = median(x, na.rm = TRUE),
        median_P2   = median(y, na.rm = TRUE),
        median_diff = median(x - y, na.rm = TRUE),
        p_value     = wt$p.value,
        statistic   = wt$statistic,
        direction   = ifelse(median(x) > median(y), "P1 > P2", "P2 > P1")
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
  P1_col <- paste0(m, "_P1")
  P2_col <- paste0(m, "_P2")
  x <- paired[[P1_col]]; y <- paired[[P2_col]]
  valid <- complete.cases(x, y); x <- x[valid]; y <- y[valid]
  if (length(x) < 3) return(NULL)
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
  tibble(
    tool_category = "ALL (pooled)", metric = m,
    n_tools = length(x),
    median_P1 = median(x), median_P2 = median(y),
    median_diff = median(x - y),
    p_value = wt$p.value, statistic = wt$statistic,
    direction = ifelse(median(x) > median(y), "P1 > P2", "P2 > P1"),
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
message("  WILCOXON SIGNED-RANK TEST: P1 (Synthetic) vs P2 V2 (Real Life)")
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
        "    %-20s  n=%-3d  P1 med=%.3f  P2 med=%.3f  diff=%.3f  W=%.0f  p=%.4f %s  [%s]",
        r$tool_category, r$n_tools, r$median_P1, r$median_P2,
        r$median_diff, r$statistic, r$p_value, r$sig_star, r$direction
      ))
    }
  }
}

message("\n", paste(rep("\u2550", 70), collapse = ""))

# ============================================================
# 8. PLOTS
# ============================================================

# Nature-mode toggle: set TRUE to suppress plot.title/plot.subtitle (titles + subtitles
# move into the figure caption text per Nature style); leave FALSE for thesis-style draft
# figures with in-figure titles. Standardised with the P1/P2 Plots_Analysis.R scripts.
NATURE_MODE <- T

theme_bench <- function(base_size = 8, nature_mode = NATURE_MODE) {
  theme_minimal(base_size = base_size) +
    theme(
      text               = element_text(family = "Helvetica"),
      plot.title         = if (nature_mode) element_blank() else element_text(face = "bold", size = base_size + 1,
                                                                              colour = "grey10", margin = margin(b = 4)),
      plot.subtitle      = if (nature_mode) element_blank() else element_text(colour = "grey40", size = base_size - 1,
                                                                             margin = margin(b = 6)),
      strip.text         = element_text(face = "bold", size = base_size - 1,
                                        margin = margin(3, 3, 3, 3)),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      panel.border       = element_rect(colour = "grey20", fill = NA, linewidth = 0.4),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line          = element_blank(),
      axis.title         = element_text(size = base_size),
      axis.text          = element_text(colour = "grey30", size = base_size - 1),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 1),
      legend.key.size    = unit(3, "mm"),
      plot.caption       = element_text(hjust = 0, size = base_size - 1,
                                        colour = "grey30",
                                        margin = margin(t = 8),
                                        lineheight = 1.1),
      plot.margin        = margin(6, 8, 16, 8)
    )
}

# --- 8a. Paired difference dot plot (kappa and F1) ---
paired_long <- paired %>%
  pivot_longer(
    cols      = c(kappa_P1, kappa_P2),
    names_to  = c("metric", "pipeline"),
    names_pattern = "(kappa)_(P1|P2)",
    values_to = "value"
  ) %>%
  mutate(
    pipeline = recode(pipeline, P1 = "Synthetic (P1)", P2 = "Real Life (P2)"),
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
    values = c("Synthetic (P1)" = "#3498DB", "Real Life (P2)" = "#E74C3C"),
    name   = "Pipeline"
  ) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Paired tool performance: Synthetic (P1) vs Real Life (P2 V2)",
    subtitle = "Each line connects the same tool across pipelines. Faceted by category and metric.",
    x = "Score", y = NULL
  ) +
  theme_bench() +
  theme(axis.text.y = element_text(size = 7),
        strip.text.y = element_text(angle = 0, hjust = 0))

# --- 8b. Category-level summary data for significance stars ---
diff_data <- paired %>%
  mutate(
    kappa_diff = kappa_P1 - kappa_P2,
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
    title    = "Per-tool score change: Synthetic (P1) \u2212 Real Life (P2 V2)",
    subtitle = "Each point = one tool. Values above zero = better on synthetic data. Stars = Wilcoxon signed-rank p-value.",
    x = NULL, y = "Score difference (P1 \u2212 P2)"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Render and save ---
plot(p_paired)
plot(p_box)

ggsave("P1vP2_wilcoxon_paired_dotplot.png",    p_paired,  width = 356, height = 406, units = "mm", dpi = 1200, path = PLOT_DIR)
ggsave("P1vP2_wilcoxon_diff_boxplot.png",      p_box,     width = 305, height = 152, units = "mm", dpi = 1200, path = PLOT_DIR)

message("\nAll plots saved to: ", PLOT_DIR)

# ============================================================
# 9. EXPORT RESULTS TABLE
# ============================================================

write.csv(all_results, file.path(PLOT_DIR, "P1vP2_wilcoxon_results_table.csv"), row.names = FALSE)
message("Results table saved to: ", file.path(PLOT_DIR, "wilcoxon_results_table.csv"))

# ============================================================
# 10. EXPORT CATEGORY AVERAGE KAPPA TABLE
# ============================================================

# Average kappa for each tool category in each P1 scenario
P1_cat_scenario <- df_P1 %>%
  dplyr::filter(!is.na(tool_category)) %>%
  dplyr::group_by(tool_category, scenario) %>%
  dplyr::summarise(kappa_mean = mean(kappa_mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(scenario = paste0("P1_", scenario)) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = kappa_mean)

# Average kappa for each tool category in each P2 scenario
P2_cat_scenario <- df_P2 %>%
  dplyr::filter(!is.na(tool_category)) %>%
  dplyr::mutate(scenario = sub("-.*", "", scenario)) %>% # Extract 'S1', 'S2', etc.
  dplyr::group_by(tool_category, scenario) %>%
  dplyr::summarise(kappa_mean = mean(kappa_mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(scenario = paste0("P2_", scenario)) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = kappa_mean)

# Average across all of P1 per tool category
P1_cat_overall <- df_P1 %>%
  dplyr::filter(!is.na(tool_category)) %>%
  dplyr::group_by(tool_category) %>%
  dplyr::summarise(P1_Overall = mean(kappa_mean, na.rm = TRUE), .groups = "drop")

# Average across all of P2 per tool category
P2_cat_overall <- df_P2 %>%
  dplyr::filter(!is.na(tool_category)) %>%
  dplyr::group_by(tool_category) %>%
  dplyr::summarise(P2_Overall = mean(kappa_mean, na.rm = TRUE), .groups = "drop")

# Merge into a single table
category_kappa_table <- P1_cat_scenario %>%
  dplyr::left_join(P2_cat_scenario, by = "tool_category") %>%
  dplyr::left_join(P1_cat_overall, by = "tool_category") %>%
  dplyr::left_join(P2_cat_overall, by = "tool_category")

write.csv(category_kappa_table, file.path(PLOT_DIR, "P1vP2_category_kappa_averages.csv"), row.names = FALSE)
message("Category average kappa table saved to: ", file.path(PLOT_DIR, "P1vP2_category_kappa_averages.csv"))
