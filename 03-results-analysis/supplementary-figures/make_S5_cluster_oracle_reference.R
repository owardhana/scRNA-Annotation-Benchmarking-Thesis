# ============================================================
# Supplementary Fig. S5 — Cluster-oracle reference heatmap
# CIPR + scmap (cluster-level) under the R+C oracle.
# Two arms: Synthetic (Phase 1 = results/P4) and Real (Phase 2 = results/P5).
#
# These two correlation-based tools are EXCLUDED from the main
# P4/P5 tool heatmaps (exclude_tools list in each Plots_Analysis.R),
# so this script re-loads the same per-scenario CSVs WITHOUT the
# exclusion and keeps only CIPR + scmap_cluster.
#
# Output: plot_S5_cluster_oracle_reference.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})

# RESULTS should point at the original (non-archived) results/ tree containing the phase-level
# per-scenario CSVs this script re-loads (see P4/P5_* below); not included in this archive.
RESULTS <- Sys.getenv("THESIS_RESULTS_DIR", "results")
OUT_DIR <- file.path(RESULTS, "_SI_generated")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

KEEP_TOOLS <- c("CIPR", "scmap_cluster")

# --- house theme (copied from the phase Plots_Analysis.R) ---
theme_bench <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(family = "Helvetica"),
      plot.title       = element_text(face = "bold", size = base_size + 2, colour = "grey10", margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey40", size = base_size, margin = margin(b = 6)),
      strip.text       = element_text(face = "bold", size = base_size, margin = margin(3,3,3,3)),
      strip.background = element_rect(fill = "grey92", colour = NA),
      panel.border     = element_rect(colour = "grey20", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      plot.caption     = element_text(hjust = 0, size = base_size - 1, colour = "grey30",
                                      margin = margin(t = 8), lineheight = 1.1),
      plot.margin      = margin(6, 8, 10, 8)
    )
}

.drop_index_col <- function(df) {
  if (ncol(df) > 0 && (names(df)[1] == "" || names(df)[1] %in% c("X", "X.1"))) df <- df[, -1, drop = FALSE]
  df
}

# Read every CSV under a scenario folder (recursing into rep dirs), keep KEEP_TOOLS.
load_arm <- function(base_dir, subfolders, scen_ids) {
  out <- list()
  for (i in seq_along(subfolders)) {
    sp <- file.path(base_dir, subfolders[i])
    if (!dir.exists(sp)) { warning("missing folder: ", sp); next }
    csvs <- list.files(sp, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
    rows <- lapply(csvs, function(f) {
      tryCatch({
        df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
        df <- .drop_index_col(df)
        if (!"tool" %in% names(df)) return(NULL)
        df <- df[df$tool %in% KEEP_TOOLS, , drop = FALSE]
        if (!nrow(df)) return(NULL)
        km <- if ("kappa_mean" %in% names(df)) df$kappa_mean else df$pooled_kappa
        pk <- if ("pooled_kappa" %in% names(df)) df$pooled_kappa else km
        am <- if ("accuracy_mean" %in% names(df)) df$accuracy_mean else NA
        fm <- if ("f1_mean" %in% names(df)) df$f1_mean else NA
        data.frame(tool = df$tool, kappa_mean = as.numeric(km), pooled_kappa = as.numeric(pk),
                   accuracy_mean = as.numeric(am), f1_mean = as.numeric(fm),
                   stringsAsFactors = FALSE)
      }, error = function(e) NULL)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) { d <- bind_rows(rows); d$scen <- scen_ids[i]; out[[subfolders[i]]] <- d }
  }
  bind_rows(out)
}

# same NA rule as the phase scripts: degenerate (0/0/0) -> NA
apply_na_rule <- function(d) {
  d %>% mutate(kappa_use = ifelse(!is.na(accuracy_mean) & !is.na(f1_mean) &
                                    kappa_mean == 0 & accuracy_mean == 0 & f1_mean == 0,
                                  NA_real_, kappa_mean),
               kappa_use = ifelse(is.na(kappa_use), pooled_kappa, kappa_use))
}

scen_ids <- paste0("S", 1:9)

synth <- load_arm(
  file.path(RESULTS, "P4 Synthetic Datasets V1"),
  paste0("S", 1:9), scen_ids) %>% apply_na_rule()

real <- load_arm(
  file.path(RESULTS, "P5 Real Life Datasets V2"),
  c("S1-Darmanis-Brain-2015_for_use","S2-Marques-Brain_for_use",
    "S3-Nowakowski-Cortex-2017_for_use","S4-Grun-Pancreas-2016_original",
    "S5-Tabula Muris-FACS-3k_for_use","S6-He-Skin-2020_for_use",
    "S7-Zhao-Immune-Fine-2020_for_use","S8-Zheng-Zhengsort-5cl-2017_original",
    "S9-MacParland-Liver-Broad_original"),
  scen_ids) %>% apply_na_rule()

summ <- function(d, arm) {
  d %>% group_by(scen, tool) %>%
    summarise(kappa = mean(kappa_use, na.rm = TRUE), .groups = "drop") %>%
    mutate(kappa = ifelse(is.nan(kappa), NA_real_, kappa), arm = arm)
}

plot_df <- bind_rows(summ(synth, "Synthetic (Phase 1)"),
                     summ(real,  "Real validation (Phase 2)")) %>%
  mutate(scen = factor(scen, levels = scen_ids),
         tool = factor(tool, levels = c("scmap_cluster", "CIPR")),
         arm  = factor(arm, levels = c("Synthetic (Phase 1)", "Real validation (Phase 2)")))

cat("\n--- S5 data (mean kappa) ---\n"); print(as.data.frame(plot_df))

p <- ggplot(plot_df, aes(x = scen, y = tool, fill = kappa)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(kappa), "–", sprintf("%.2f", kappa))),
            size = 3.4, colour = "white", fontface = "bold") +
  facet_wrap(~ arm, ncol = 1) +
  scale_fill_gradient2(low = "#c0392b", mid = "#f39c12", high = "#27ae60",
                       midpoint = 0.5, limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1), name = "Mean Cohen's κ",
                       na.value = "grey85",
                       guide = guide_colorbar(barwidth = unit(9, "cm"), barheight = unit(0.4, "cm"),
                                              title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(labels = function(x) paste0(x, " †")) +
  labs(title = "Cluster-oracle reference: CIPR and cluster-level scmap",
       subtitle = "Mean Cohen's κ under the R+C oracle. Two correlation-based tools excluded from the Phase 1–3 cell-level ranking.",
       x = "Scenario / dataset ID (S1–S9, paired across arms)", y = "Tool",
       caption = "† Cluster-dependent tool: supplied ground-truth labels as the cluster argument (Phases 1–3).\nSynthetic κ = mean across three replicate seeds; real κ = single evaluation. – = no valid prediction.") +
  theme_bench() +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))

out <- file.path(OUT_DIR, "plot_S5_cluster_oracle_reference.png")
ggsave(out, p, width = 240, height = 130, units = "mm", dpi = 600)
cat("\nSaved:", out, "\n")
