# ============================================================
# scRNA-seq Benchmark Pipeline & Statistical Analysis
# Phase 3: Cross-Platform Batch Effect Analysis (P9)
# Merged Script: Loads data, computes statistics, and generates plots
# Run top-to-bottom in RStudio
# ============================================================
#
# INTERPRETATION NOTES:
#   - PX-1 and PX-2 near-zero kappa is expected and is a scientific
#     finding, not a pipeline error — see methods for the 5 compounding
#     reasons (UMI vs read counts, gene detection rate, T2D donors,
#     donor batch effects, asymmetric class imbalance)
#   - PBMC-B CD14+ monocyte F1 is confounded by inDrops enrichment
#     artifact — do not interpret low F1 here as tool failure
#   - PBMC-C has 6 cell types only (CD16+ mono and NK excluded due to
#     structural zeros in Seq-Well)
#   - Tools with U% = NA have no rejection option — treat as
#     structurally missing, not zero
# ============================================================

# ============================================================
# 1. SETUP & LIBRARIES
# ============================================================
# install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "patchwork",
#                    "scales", "forcats", "PMCMRplus", "ComplexHeatmap",
#                    "circlize"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(scales)
library(forcats)
library(ComplexHeatmap)
library(circlize)
library(PMCMRplus)   # Nemenyi post-hoc
library(grid)
library(extrafont)   # Arial font support for cairo_pdf
library(Cairo)       # cairo_pdf device

# Cluster-dependent tool decoration (dagger suffix on tool labels)
# Defines CLUSTER_DEPENDENT_TOOLS, mark_cluster_dep(), CLUSTER_DEP_CAPTION
source("../_cluster_dependent.R")

# ============================================================
# 2. CONSTANTS & CONFIGURATIONS
# ============================================================

# Directories
BASE_DIR   <- "./data"
PLOT_DIR   <- "./Plots"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

# Category colour palette — consistent with P4/P5/P6/P7/P8 scripts
CATEGORY_COLOURS <- c(
  "Marker-Based"      = "#E69F00",
  "Correlation-Based" = "#56B4E9",
  "Classic ML"        = "#009E73",
  "Deep Learning"     = "#CC79A7",
  "Semi-Supervised"   = "#0072B2"
)

# Trial labels for axis display
TRIAL_LABELS <- c(
  "CB-1"   = "CB-1: 10xv2 \u2192 CEL-seq2",
  "CB-2"   = "CB-2: 10xv2 \u2192 Drop-seq",
  "CB-3"   = "CB-3: CEL-seq2 \u2192 10xv2",
  "CB-4"   = "CB-4: Drop-seq \u2192 10xv2",
  "PX-1"   = "PX-1: Baron \u2192 Segerstolpe\u2020",
  "PX-2"   = "PX-2: Segerstolpe \u2192 Baron\u2020",
  "PBMC-A" = "PBMC-A: 10xv2 \u2192 Drop-seq",
  "PBMC-B" = "PBMC-B: 10xv2 \u2192 InDrops*",
  "PBMC-C" = "PBMC-C: 10xv2 \u2192 Seq-Well"
)

# Trial-to-block mapping
TRIAL_BLOCKS <- c(
  "CB-1" = "CellBench", "CB-2" = "CellBench",
  "CB-3" = "CellBench", "CB-4" = "CellBench",
  "PX-1" = "Pancreas",  "PX-2" = "Pancreas",
  "PBMC-A" = "PBMCbench", "PBMC-B" = "PBMCbench",
  "PBMC-C" = "PBMCbench"
)

# Platform shift classification per trial
PLATFORM_SHIFTS <- c(
  "CB-1" = "droplet and plate",    "CB-2" = "within droplet",
  "CB-3" = "droplet and plate",    "CB-4" = "within droplet",
  "PX-1" = "droplet and plate",    "PX-2" = "droplet and plate",
  "PBMC-A" = "within droplet",    "PBMC-B" = "within droplet",
  "PBMC-C" = "droplet and nanowell"
)

# Tool Category Lookup Vectors
marker_tools <- c(
  "clustifyr_hyper",
  "clustifyr_jaccard",
  "scCATCH",
  "SCINA",
  "ScInfeR",
  "SCSA",
  "scSorter",
  "scType"
)

correlation_tools <- c(
  "CIPR",
  "scibetR",
  "scmap_cell",
  "scmap_cluster",
  "Seurat_Transfer_CCA",
  "Seurat_Transfer_PCA",
  "Seurat_Transfer_RPCA",
  "SingleR"
)

classicML_tools <- c(
  "CaSTLe",
  "CellTypist",
  "CHETAH",
  "scAnnotate",
  "scAnnotatR",
  "scClassify",
  "scID",
  "scPred",
  "scPred_adaboost",
  "scPred_avNNet",
  "scPred_bayesglm",
  "scPred_earth",
  "scPred_glm",
  "scPred_glmboost",
  "scPred_knn",
  "scPred_lda",
  "scPred_mlp",
  "scPred_multinom",
  "scPred_nb",
  "scPred_nnet",
  "scPred_regLogistic",
  "scPred_rf",
  "scPred_svmLinear",
  "scPred_svmPoly",
  "scPred_xgbTree",
  "singleCellNet"
)

DL_tools <- c(
  "ACTINN",
  "Cell_BLAST",
  "CIForm",
  "NeuCA",
  "scBalance",
  "scDeepSort",
  "scHash",
  "scLearn",
  "scMMT",
  "TOSICA"
)

semi_supervised_tools <- c(
  "CALLR",
  "CAMLU",
  "MARS",
  "SCANVI",
  "scArches",
  "scGAD",
  "scnym",
  "scSemiCluster",
  "scSemiGAN"
)

# Returns a named character vector mapping every tool name to its category label.
.build_category_map <- function() {
  c(
    setNames(rep("Marker-Based",      length(marker_tools)),           marker_tools),
    setNames(rep("Correlation-Based", length(correlation_tools)),      correlation_tools),
    setNames(rep("Classic ML",        length(classicML_tools)),        classicML_tools),
    setNames(rep("Deep Learning",     length(DL_tools)),               DL_tools),
    setNames(rep("Semi-Supervised",   length(semi_supervised_tools)),  semi_supervised_tools)
  )
}

# Shared Theme Component — consistent with P4/P5/P6/P7/P8 scripts
theme_bench <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      text               = element_text(family = "Helvetica"),
      plot.title         = element_text(face = "bold", size = base_size + 1,
                                        colour = "grey10", margin = margin(b = 4)),
      plot.subtitle      = element_text(colour = "grey40", size = base_size - 1,
                                        margin = margin(b = 6)),
      strip.text         = element_text(face = "bold", size = base_size - 1,
                                        margin = margin(3, 3, 3, 3)),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      panel.border       = element_rect(colour = "grey20", fill = NA, linewidth = 0.4),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
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

kappa_axis   <- scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1))
kappa_axis_x <- scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1))

# ============================================================
# 3. HELPER FUNCTIONS & DATA LOADING
# ============================================================
load_benchmark_results <- function(base_dir, exclude_tools = NULL) {
  if (!dir.exists(base_dir)) stop("base_dir does not exist: ", base_dir)
  
  csv_files <- list.files(base_dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  if (length(csv_files) == 0) stop("No CSV files found in ", base_dir)
  
  cat_map <- .build_category_map()
  all_dfs <- list()
  
  for (f in csv_files) {
    fname <- basename(f)
    # Extract trial name like CB-1, PBMC-A, PX-2
    trial_match <- regexpr("(CB-[1-4]|PX-[1-2]|PBMC-[A-C])", fname)
    if (trial_match == -1) next
    trial <- regmatches(fname, trial_match)
    
    tryCatch({
      df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      
      # Determine if tool column is 1 or 2
      if ("tool" %in% colnames(df)) {
        tool_col <- "tool"
      } else if (colnames(df)[2] == "tool" || colnames(df)[1] == "V2") {
        # Edge case handling for unnamed index columns
        if(colnames(df)[2] == "tool") df <- df[, -1, drop=FALSE]
      }
      
      # Ensure consistent tool naming if tool column exists
      if ("tool" %in% colnames(df)) {
        df[["trial"]] <- trial
        all_dfs[[length(all_dfs) + 1]] <- df
      }
    }, error = function(e) { NULL })
  }
  
  combined <- dplyr::bind_rows(all_dfs)
  
  # Average over multiple runs per tool per trial if present
  combined <- combined %>%
    group_by(trial, tool) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  if (!is.null(exclude_tools) && length(exclude_tools) > 0) {
    combined <- combined[!combined[["tool"]] %in% exclude_tools, , drop = FALSE]
  }
  
  combined[["tool_category"]] <- cat_map[combined[["tool"]]]
  
  combined <- combined %>% filter(!is.na(tool_category))
  combined[["tool_category"]] <- factor(combined[["tool_category"]], levels = names(CATEGORY_COLOURS))
  combined[["block"]] <- factor(TRIAL_BLOCKS[combined$trial], levels = c("CellBench", "Pancreas", "PBMCbench"))
  combined[["platform_shift"]] <- factor(PLATFORM_SHIFTS[combined$trial], levels = c("within droplet", "droplet and nanowell", "droplet and plate"))
  
  message("Loaded ", nrow(combined), " tool-trial pairs")
  combined
}

# Load data
df <- load_benchmark_results(
  base_dir = BASE_DIR,
  exclude_tools = c('scInfeR', 'scPred_adaboost', "scGAD", "scDeepSort", "CAMLU", "scAnnotatR")
)

df <- df %>%
  mutate(kappa_mean = ifelse(kappa_mean == 0 & accuracy_mean == 0 & f1_mean == 0, NA_real_, kappa_mean))

# Expand to full grid (all tools x all trials)
all_tools <- unique(df$tool)
all_trials <- names(TRIAL_BLOCKS)
df_full <- expand_grid(tool = all_tools, trial = all_trials) %>%
  left_join(df, by = c("tool", "trial")) %>%
  mutate(
    tool_category = .build_category_map()[tool],
    tool_category = factor(tool_category, levels = names(CATEGORY_COLOURS)),
    block = factor(TRIAL_BLOCKS[trial], levels = c("CellBench", "Pancreas", "PBMCbench")),
    platform_shift = factor(PLATFORM_SHIFTS[trial], levels = c("within droplet", "droplet and nanowell", "droplet and plate"))
  )

# ============================================================
# 4. SECTION 1: DATA LOADING AND VALIDATION
# ============================================================

# Tool failure summary table
# A tool fails if it has NA for kappa_mean
failure_summary <- df_full %>%
  mutate(status = case_when(
    is.na(kappa_mean) ~ "FAIL/NA",
    TRUE ~ "PASS"
  )) %>%
  dplyr::select(tool, trial, status) %>%
  pivot_wider(names_from = trial, values_from = status)


# ============================================================
# 5. SECTION 2: OVERVIEW HEATMAPS
# ============================================================

make_heatmap <- function(metric_col, title, is_unassigned=FALSE) {

  legend_name <- paste("Mean", title)

  plot_data <- df_full %>%
    dplyr::select(tool, tool_category, trial, block, platform_shift, !!sym(metric_col)) %>%
    group_by(tool) %>%
    mutate(overall_metric = mean(!!sym(metric_col), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      tool = factor(tool, levels = rev(names(.build_category_map()))),
      trial = factor(trial, levels = all_trials)
    )

  fill_guide <- guide_colorbar(barwidth = unit(10, "cm"),
                               barheight = unit(0.4, "cm"),
                               title.position = "top",
                               title.hjust = 0.5)

  if (is_unassigned) {
    fill_scale <- scale_fill_gradient2(
      low = "#27ae60", mid = "#f39c12", high = "#c0392b",
      midpoint = 0.5, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = legend_name, na.value = "white",
      guide = fill_guide
    )
  } else {
    fill_scale <- scale_fill_gradient2(
      low = "#c0392b", mid = "#f39c12", high = "#27ae60",
      midpoint = 0.5, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = legend_name, na.value = "white",
      guide = fill_guide
    )
  }

  annot_data <- plot_data %>% 
    distinct(trial, block, platform_shift) %>%
    mutate(dummy = "")

  p_annot <- ggplot(annot_data, aes(x = trial, y = "Involved Platforms", fill = platform_shift)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(dummy ~ block, scales = "free", space = "free") +
    scale_fill_manual(
      values = c("within droplet" = "#d5d8dc", "droplet and nanowell" = "#85929e", "droplet and plate" = "#2e4053"),
      name = "Involved Platforms"
    ) +
    theme_void() +
    theme(
      strip.text.x = element_blank(),
      strip.text.y = element_text(colour = "transparent", size = 8),
      axis.text.y = element_text(size = 9, face = "bold", margin = margin(r = 5), hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 10, b = -5, l = 10, r = 10)
    )

  p_main <- ggplot(plot_data, aes(x = trial, y = tool, fill = !!sym(metric_col))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(tool_category ~ block, scales = "free", space = "free") +
    fill_scale +
    scale_y_discrete(labels = mark_cluster_dep) +
    labs(
      title    = paste("Per-tool", title, "across trials (Cross-Platform Validation)"),
      subtitle = "Tools sorted alphabetically within each category. Faceted by algorithm family and trial block.",
      x        = "Trial",
      y        = "Tool",
      caption  = CLUSTER_DEP_CAPTION
    ) +
    theme_bench() +
    theme(
      axis.text.y      = element_text(size = 8),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      strip.text.y     = element_text(angle = 0, hjust = 0),
      strip.text.x     = element_text(face = "bold", size = 10),
      panel.grid.major = element_blank(),
      plot.margin      = margin(t = 5, b = 16, l = 12, r = 12),
      plot.caption     = element_text(hjust = 0, size = 8, colour = "grey30",
                                      margin = margin(t = 8))
    )

  p_combined <- p_annot / p_main + plot_layout(heights = c(1, 25), guides = "collect")
  return(p_combined)
}

# Plot 1: Kappa Heatmap
p1 <- make_heatmap("kappa_mean", "Cohen's \u03ba")

# Plot 2: Macro F1 Heatmap
p2 <- make_heatmap("f1_mean", "Macro F1")

# Plot 3: Unassigned Rate Heatmap
p3 <- make_heatmap("unassigned_mean", "Unassigned %", is_unassigned=TRUE)

# ============================================================
# 6. SECTION 3: WITHIN-BLOCK COMPARISONS
# ============================================================

# Plot 8: CellBench Direction Effect
cb_data <- df_full %>%
  filter(block == "CellBench", !is.na(kappa_mean)) %>%
  mutate(
    direction = case_when(
      trial %in% c("CB-1", "CB-2") ~ "10X -> other",
      trial %in% c("CB-3", "CB-4") ~ "other -> 10X"
    ),
    pair = case_when(
      trial %in% c("CB-1", "CB-3") ~ "CEL-seq2",
      trial %in% c("CB-2", "CB-4") ~ "Drop-seq"
    )
  )

# Highlight tools where direction reversal causes >0.1 drop
cb_diff <- cb_data %>%
  dplyr::select(tool, tool_category, pair, direction, kappa_mean) %>%
  pivot_wider(names_from = direction, values_from = kappa_mean) %>%
  mutate(diff = `10X -> other` - `other -> 10X`,
         highlight = abs(diff) > 0.1)

cb_plot_data <- cb_data %>% left_join(cb_diff %>% dplyr::select(tool, pair, highlight), by=c("tool", "pair"))

p4 <- ggplot(cb_plot_data, aes(x = direction, y = kappa_mean, group = tool, colour = tool_category)) +
  geom_line(aes(alpha = highlight, linewidth = highlight)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~pair) +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
  scale_linewidth_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  kappa_axis +
  labs(
    title    = "CellBench transfer direction asymmetry (Cross-Platform Validation)",
    subtitle = "Lines connect same tool across directions. Bold lines = >0.1 \u03ba difference.",
    x        = "Transfer direction",
    y        = "Cohen\u2019s \u03ba"
  ) +
  theme_bench()


# Plot 9: PBMC Platform Gradient
pbmc_data <- df_full %>% filter(block == "PBMCbench", !is.na(kappa_mean))
p5 <- ggplot(pbmc_data, aes(x = trial, y = kappa_mean, fill = tool_category, colour = tool_category)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, coef = 1.5) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.7) +
  scale_fill_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  kappa_axis +
  annotate("text", x = 2, y = 1.02, label = "* PBMC-B: CD14+ mono enrichment artifact in InDrops",
           fontface = "italic", size = 3, colour = "grey40") +
  labs(
    title    = "PBMCbench platform shift gradient (Cross-Platform Validation)",
    subtitle = "Each point = one tool. Boxes span IQR. PBMC-C excludes CD16+ mono and NK (structural zeros).",
    x        = "Trial",
    y        = "Cohen\u2019s \u03ba"
  ) +
  theme_bench()


# Plot 10: Pancreas Pair
# NOTE: PX-1 and PX-2 near-zero kappa is expected — UMI vs read counts,
# gene detection rate, T2D donors, donor batch effects, asymmetric class imbalance
px_data <- df_full %>% filter(block == "Pancreas", !is.na(kappa_mean))
p6 <- ggplot(px_data, aes(x = trial, y = kappa_mean, colour = tool_category)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "#c0392b", linewidth = 0.6) +
  geom_text_repel(data = filter(px_data, kappa_mean > 0.4), aes(label = mark_cluster_dep(tool)),
                  size = 2.8, max.overlaps = 20, segment.colour = "grey60", segment.size = 0.3) +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  kappa_axis +
  labs(
    title    = "Pancreas batch effect severity (Cross-Platform Validation)",
    # NOTE: the pre-existing "\u2020" prefix on "Hardest transfer pair" was dropped
    # to avoid collision with the cluster-dependent dagger convention.
    subtitle = "Hardest transfer pair. Dashed line = near-failure threshold (\u03ba = 0.2). Labels shown for \u03ba > 0.4.",
    x        = "Trial",
    y        = "Cohen\u2019s \u03ba",
    caption  = CLUSTER_DEP_CAPTION
  ) +
  theme_bench() + CLUSTER_DEP_THEME


# ============================================================
# 7. SECTION 4: CROSS-BLOCK COMPARISON
# ============================================================

# Plot 11: Batch Robustness Gradient Plot
robustness_data <- df_full %>%
  filter(!is.na(kappa_mean)) %>%
  group_by(tool, tool_category, platform_shift) %>%
  summarise(mean_kappa = mean(kappa_mean), .groups = "drop")

# Highlight flattest/steepest
slope_calcs <- robustness_data %>%
  pivot_wider(names_from = platform_shift, values_from = mean_kappa) %>%
  mutate(drop = `within droplet` - `droplet and plate`) %>%
  filter(!is.na(drop))

steepest <- slope_calcs %>% top_n(3, drop) %>% pull(tool)
flattest <- slope_calcs %>% top_n(-3, drop) %>% pull(tool)

robustness_plot_data <- robustness_data %>%
  mutate(highlight = case_when(
    tool %in% steepest ~ "Steepest decline",
    tool %in% flattest ~ "Most robust",
    TRUE ~ "Normal"
  ))

p7 <- ggplot(robustness_plot_data, aes(x = platform_shift, y = mean_kappa, group = tool, colour = tool_category)) +
  geom_line(aes(linewidth = highlight, alpha = highlight)) +
  geom_point(aes(size = highlight)) +
  scale_colour_manual(values = CATEGORY_COLOURS) +
  scale_linewidth_manual(values = c("Normal" = 0.5, "Most robust" = 1.5, "Steepest decline" = 1.5)) +
  scale_alpha_manual(values = c("Normal" = 0.3, "Most robust" = 1, "Steepest decline" = 1)) +
  scale_size_manual(values = c("Normal" = 1, "Most robust" = 3, "Steepest decline" = 3)) +
  geom_text_repel(data = filter(robustness_plot_data, platform_shift == "droplet and plate" & highlight != "Normal"),
                  aes(label = mark_cluster_dep(tool)), nudge_x = 0.2, direction = "y", hjust = 0, size = 2.8,
                  segment.colour = "grey60", segment.size = 0.3) +
  labs(
    title    = "Batch robustness gradient \u2014 slope plot (Cross-Platform Validation)",
    subtitle = "Mean \u03ba per platform shift severity. Bold lines = top 3 steepest decline / most robust tools.",
    x        = "Platform shift severity (easy \u2192 hard)",
    y        = "Mean Cohen\u2019s \u03ba",
    caption  = CLUSTER_DEP_CAPTION
  ) +
  theme_bench() + CLUSTER_DEP_THEME


# Plot 12: Biology Difficulty vs Platform Shift
scatter2d_data <- df_full %>%
  filter(!is.na(kappa_mean)) %>%
  group_by(tool, tool_category, block) %>%
  summarise(mean_kappa = mean(kappa_mean), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = mean_kappa) %>%
  filter(!is.na(CellBench), !is.na(Pancreas), !is.na(PBMCbench))

p8 <- ggplot(scatter2d_data, aes(x = CellBench, y = Pancreas, colour = tool_category, size = PBMCbench)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = mark_cluster_dep(tool)), size = 2.8, max.overlaps = 20,
                  segment.colour = "grey60", segment.size = 0.3) +
  annotate("text", x = 0.92, y = 0.95, label = "Robust",
           hjust = 1, vjust = 0, fontface = "italic", colour = "grey50", size = 3.5) +
  annotate("text", x = 0.08, y = 0.05, label = "Fails everywhere",
           hjust = 0, vjust = 1, fontface = "italic", colour = "grey50", size = 3.5) +
  annotate("text", x = 0.08, y = 0.95, label = "Biology-dependent",
           hjust = 0, vjust = 0, fontface = "italic", colour = "grey50", size = 3.5) +
  annotate("text", x = 0.92, y = 0.05, label = "Platform-sensitive",
           hjust = 1, vjust = 1, fontface = "italic", colour = "grey50", size = 3.5) +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_size_continuous(name = "Mean \u03ba\nPBMCbench", range = c(2, 7)) +
  kappa_axis + kappa_axis_x +
  labs(
    title    = "Biology difficulty vs platform shift (Cross-Platform Validation)",
    subtitle = "Dashed line = equal performance. Point size \u221d mean \u03ba on PBMCbench (medium difficulty).",
    x        = "Mean Cohen\u2019s \u03ba on CellBench (easy biology)",
    y        = "Mean Cohen\u2019s \u03ba on Pancreas (hard biology + hard platform)",
    caption  = CLUSTER_DEP_CAPTION
  ) +
  theme_bench() + CLUSTER_DEP_THEME


# ============================================================
# 8. SECTION 5: CATEGORY-LEVEL ANALYSIS
# ============================================================

# Plot 13: Category Boxplots 3x3 Grid
p9 <- df_full %>% filter(!is.na(kappa_mean)) %>%
  ggplot(aes(x = tool_category, y = kappa_mean, fill = tool_category, colour = tool_category)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, coef = 1.5) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.7) +
  facet_wrap(~trial, ncol = 3) +
  scale_fill_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  kappa_axis +
  labs(
    title    = "Cohen\u2019s \u03ba by algorithm category and trial (Cross-Platform Validation)",
    subtitle = "Each point = one tool. Boxes span IQR. Coloured by algorithm category.",
    x        = "Tool category",
    y        = "Mean Cohen\u2019s \u03ba"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Table 14: Category Ranking Table
cat_ranking <- df_full %>%
  filter(!is.na(kappa_mean)) %>%
  group_by(block, trial, tool_category) %>%
  summarise(
    med_kappa = median(kappa_mean),
    iqr_kappa = IQR(kappa_mean),
    val_str = sprintf("%.3f \u00b1 %.3f", med_kappa, iqr_kappa),
    .groups = "drop"
  ) %>%
  dplyr::select(block, trial, tool_category, val_str) %>%
  pivot_wider(names_from = trial, values_from = val_str)


# Plot 15: Category Robustness Violin (Proxy)
# Proxy: max kappa CellBench - kappa Pancreas
robust_proxy <- df_full %>% filter(!is.na(kappa_mean)) %>% group_by(tool, tool_category) %>%
  summarise(
    max_cb = max(kappa_mean[block == "CellBench"], na.rm = TRUE),
    mean_px = mean(kappa_mean[block == "Pancreas"], na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(is.finite(max_cb), is.finite(mean_px)) %>%
  mutate(kappa_drop = max_cb - mean_px)

p10 <- ggplot(robust_proxy, aes(x = tool_category, y = kappa_drop, fill = tool_category, colour = tool_category)) +
  geom_violin(alpha = 0.15, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Category robustness under batch effects (Cross-Platform Validation)",
    subtitle = "Proxy: max \u03ba on CellBench \u2212 mean \u03ba on Pancreas. Higher = more degradation.",
    x        = "Tool category",
    y        = "\u03ba drop"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# ============================================================
# 9. SECTION 7: RUNTIME AND EFFICIENCY
# ============================================================

# Plot 18: Runtime vs Kappa scatter
runtime_data <- df_full %>% filter(!is.na(kappa_mean), !is.na(runtime_mean), runtime_mean > 0) %>%
  group_by(tool, tool_category, block) %>%
  summarise(mean_kappa = mean(kappa_mean), mean_rt = mean(runtime_mean), .groups = "drop") %>%
  mutate(log_rt = log10(mean_rt))

# Simple Pareto frontier function per block
get_pareto <- function(df) {
  df <- df[order(df$mean_rt, -df$mean_kappa), ]
  pareto <- df[0, ]
  best_kappa <- -1
  for(i in 1:nrow(df)) {
    if(df$mean_kappa[i] > best_kappa) {
      pareto <- rbind(pareto, df[i, ])
      best_kappa <- df$mean_kappa[i]
    }
  }
  pareto
}

pareto_fronts <- runtime_data %>% group_by(block) %>% group_modify(~get_pareto(.x)) %>% ungroup()

p11 <- ggplot(runtime_data, aes(x = mean_rt, y = mean_kappa, colour = tool_category)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_line(data = pareto_fronts, aes(x = mean_rt, y = mean_kappa),
            colour = "grey50", linetype = "dashed", inherit.aes = FALSE) +
  geom_text_repel(data = pareto_fronts, aes(label = mark_cluster_dep(tool)), size = 2.8,
                  show.legend = FALSE, segment.colour = "grey60", segment.size = 0.3) +
  facet_wrap(~block, ncol = 3, scales = "free_x") +
  scale_x_log10(labels = label_number(suffix = "s", big.mark = ",")) +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  kappa_axis +
  labs(
    title    = "Performance\u2013runtime tradeoff (Cross-Platform Validation)",
    subtitle = "Each point = one tool. Dashed line = Pareto frontier. Labels = Pareto-optimal tools.",
    x        = "Mean runtime (seconds, log\u2081\u2080 scale)",
    y        = "Mean Cohen\u2019s \u03ba",
    caption  = CLUSTER_DEP_CAPTION
  ) +
  theme_bench() + CLUSTER_DEP_THEME


# ============================================================
# 10. SECTION 8: STATISTICAL TESTS
# ============================================================

# 19. Friedman Test per Block + CD diagram
make_block_cd <- function(block_id) {
  d <- df_full %>% filter(block == block_id, !is.na(kappa_mean)) %>% dplyr::select(trial, tool, kappa_mean)
  mat_wide <- d %>% pivot_wider(names_from = tool, values_from = kappa_mean)
  mat <- as.matrix(mat_wide[, -1])
  rownames(mat) <- mat_wide$trial
  
  # keep complete
  complete_s <- colnames(mat)[colSums(is.na(mat)) == 0]
  if(length(complete_s) < 2 || nrow(mat) < 2) return(NULL)
  mat <- mat[, complete_s, drop=FALSE]
  
  fr_res <- friedman.test(mat)
  k_s <- ncol(mat)
  N_s <- nrow(mat)
  
  # Simplified CD
  q_alpha <- qtukey(1 - 0.05, k_s, Inf) / sqrt(2)
  CD_s <- q_alpha * sqrt(k_s * (k_s + 1) / (6 * N_s))
  
  rank_mat <- t(apply(mat, 1, function(row) rank(-row, ties.method="average")))
  mean_ranks <- colMeans(rank_mat)
  
  mr_df <- tibble(tool = names(mean_ranks), mean_rank = mean_ranks) %>% arrange(mean_rank) %>%
    left_join(df_full %>% distinct(tool, tool_category), by="tool")
  
  p_val <- fr_res$p.value
  fr_label <- sprintf("Friedman p = %.3e", p_val)
  
  p <- ggplot(mr_df, aes(x = mean_rank, y = 0)) +
    annotate("segment", x = 1, xend = 1 + CD_s, y = 0.25, yend = 0.25,
             colour = "black", linewidth = 1.5, lineend = "round") +
    annotate("text", x = 1 + CD_s/2, y = 0.29,
             label = sprintf("CD = %.2f", CD_s), size = 2.8, fontface = "bold", vjust = 0) +
    geom_point(aes(colour = tool_category), size = 3, alpha = 0.9) +
    geom_text_repel(aes(label = mark_cluster_dep(tool), colour = tool_category),
                    nudge_y = 0.14, direction = "x",
                    segment.size = 0.3, segment.colour = "grey60",
                    size = 2.2, max.overlaps = 40) +
    scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
    scale_x_continuous(name = "Mean rank (1 = best)",
                       breaks = seq(1, k_s, by = max(1, round(k_s / 10)))) +
    coord_cartesian(ylim = c(-0.28, 0.42)) +
    labs(title = paste(block_id, "\u2014 critical difference diagram"),
         subtitle = fr_label, y = NULL,
         caption = CLUSTER_DEP_CAPTION) +
    theme_bench() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    CLUSTER_DEP_THEME
  
  return(p)
}

cd_plots <- purrr::compact(lapply(unique(df_full$block), make_block_cd))
if(length(cd_plots) > 0) {
  p12 <- wrap_plots(cd_plots, ncol = 1)
}

# 20. Direction Asymmetry Test
stats_results <- list()
cb1_3 <- cb_data %>% filter(pair == "CEL-seq2") %>% dplyr::select(tool, direction, kappa_mean) %>% pivot_wider(names_from=direction, values_from=kappa_mean) %>% filter(!is.na(`10X -> other`), !is.na(`other -> 10X`))
res1 <- wilcox.test(cb1_3$`10X -> other`, cb1_3$`other -> 10X`, paired=TRUE)
stats_results[[1]] <- data.frame(Test="Direction Asymmetry CEL-seq2", p_value=res1$p.value, stat=res1$statistic)

cb2_4 <- cb_data %>% filter(pair == "Drop-seq") %>% dplyr::select(tool, direction, kappa_mean) %>% pivot_wider(names_from=direction, values_from=kappa_mean) %>% filter(!is.na(`10X -> other`), !is.na(`other -> 10X`))
res2 <- wilcox.test(cb2_4$`10X -> other`, cb2_4$`other -> 10X`, paired=TRUE)
stats_results[[2]] <- data.frame(Test="Direction Asymmetry Drop-seq", p_value=res2$p.value, stat=res2$statistic)

# 21. Platform Shift ANOVA
anova_data <- df_full %>% filter(block != "Pancreas", !is.na(kappa_mean))
for(cat in names(CATEGORY_COLOURS)) {
  d <- anova_data %>% filter(tool_category == cat)
  if(nrow(d) > 3) {
    fit <- aov(kappa_mean ~ platform_shift, data=d)
    s <- summary(fit)
    pval <- s[[1]][["Pr(>F)"]][1]
    fval <- s[[1]][["F value"]][1]
    stats_results[[length(stats_results)+1]] <- data.frame(Test=paste("ANOVA", cat), p_value=pval, stat=fval)
  }
}


# 22. Top Tool Identification
top_tools <- df_full %>% filter(!is.na(kappa_mean)) %>% group_by(trial) %>%
  arrange(desc(kappa_mean)) %>%
  slice_head(n=3) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(trial, rank, tool, kappa_mean, tool_category)


# 23. Category Scatter Barplot across Platform Shifts
shift_bar_summary <- df_full %>%
  filter(!is.na(kappa_mean)) %>%
  group_by(tool_category, platform_shift) %>%
  summarise(
    med_kappa = median(kappa_mean),
    q25 = quantile(kappa_mean, 0.25),
    q75 = quantile(kappa_mean, 0.75),
    .groups = "drop"
  )

p13 <- ggplot() +
  geom_col(data = shift_bar_summary, aes(x = platform_shift, y = med_kappa, fill = tool_category),
           position = position_dodge(width = 0.8), alpha = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(data = shift_bar_summary, aes(x = platform_shift, ymin = q25, ymax = q75, group = tool_category),
                position = position_dodge(width = 0.8), width = 0.25, colour = "black", linewidth = 0.5) +
  geom_point(data = filter(df_full, !is.na(kappa_mean)), 
             aes(x = platform_shift, y = kappa_mean, colour = tool_category),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  scale_colour_manual(values = CATEGORY_COLOURS, name = "Tool Category") +
  kappa_axis +
  labs(
    title    = "Cohen\u2019s \u03ba by tool category and platform shift",
    subtitle = "Bars = median, Error bars = IQR. Points represent individual tools.",
    x        = "Platform Shift",
    y        = "Cohen\u2019s \u03ba"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# ============================================================
# 11. SAVE OUTPUTS
# ============================================================

message("\nGenerating and saving plots to: ", PLOT_DIR)

write.csv(failure_summary, file.path(PLOT_DIR, "failure_summary.csv"), row.names = FALSE)
ggsave("Plot1_heatmap_kappa.png",            p1,  width = 279, height = 381, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot2_heatmap_f1.png",               p2,  width = 279, height = 381, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot3_heatmap_unassigned.png",       p3,  width = 279, height = 381, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot4_cellbench_direction.png",      p4,  width = 254, height = 127, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot5_pbmc_gradient.png",            p5,  width = 254, height = 127, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot6_pancreas_pair.png",            p6,  width = 203, height = 127, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot7_batch_robustness_slope.png",   p7,  width = 254, height = 152, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot8_biology_vs_platform.png",      p8,  width = 254, height = 203, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot9_category_boxplots.png",        p9,  width = 330, height = 254, units = "mm", dpi = 600, path = PLOT_DIR)
write.csv(cat_ranking, file.path(PLOT_DIR, "category_table.csv"), row.names = FALSE)
ggsave("Plot10_category_robustness.png",     p10, width = 229, height = 127, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot11_runtime_scatter.png",         p11, width = 356, height = 127, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot12_critical_difference.png",     p12, width = 356, height = 305, units = "mm", dpi = 600, path = PLOT_DIR)
ggsave("Plot13_category_shift_barplot.png",  p13, width = 254, height = 152, units = "mm", dpi = 600, path = PLOT_DIR)
write.csv(bind_rows(stats_results), file.path(PLOT_DIR, "statistical_summary.csv"), row.names = FALSE)
write.csv(top_tools, file.path(PLOT_DIR, "top_tools_per_trial.csv"), row.names = FALSE)

message("\nAll outputs saved successfully.")
