# ============================================================
# Supplementary Fig. S7 — Per-tool runtime box-plot, LLM family (Phase 4).
# Data: results/P7b LLM - Real Internet/with_score/runtime_dataset.tsv
#       (columns: dataset, tool, runtime ; one row per tool x dataset)
# Matches the P7a marker-family runtime box-plot (plot13) styling.
# Output: plot_S7_llm_runtime_boxplot.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})

# RESULTS should point at the original (non-archived) results/ tree containing the runtime TSV
# this script re-loads (see the tsv path below); not included in this archive.
RESULTS <- Sys.getenv("THESIS_RESULTS_DIR", "results")
OUT_DIR <- file.path(RESULTS, "_SI_generated")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

tsv <- file.path(RESULTS, "P7b LLM - Real Internet", "with_score", "runtime_dataset.tsv")
rt  <- read.delim(tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
names(rt) <- trimws(names(rt))
rt <- rt %>% mutate(runtime = as.numeric(runtime)) %>% filter(!is.na(runtime), runtime > 0)

theme_bench <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(family = "Helvetica"),
      plot.title       = element_text(face = "bold", size = base_size + 2, colour = "grey10", margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey40", size = base_size, margin = margin(b = 6)),
      panel.border     = element_rect(colour = "grey20", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text        = element_text(colour = "grey30", size = base_size),
      plot.caption     = element_text(hjust = 0, size = base_size - 1, colour = "grey30", margin = margin(t = 8)),
      plot.margin      = margin(6, 10, 10, 10),
      legend.position  = "none"
    )
}

# order tools by median runtime (fast -> slow): GPTCelltype, mLLMCelltype, CASSIA
ord <- rt %>% group_by(tool) %>% summarise(m = median(runtime), .groups = "drop") %>% arrange(m)
rt$tool <- factor(rt$tool, levels = ord$tool)

TOOL_COLOURS <- c(GPTCelltype = "#3498DB", mLLMCelltype = "#9B59B6", CASSIA = "#E74C3C")

cat("\n--- S7 runtime summary (s) ---\n")
print(rt %>% group_by(tool) %>%
        summarise(n = n(), min = round(min(runtime)), median = round(median(runtime)),
                  max = round(max(runtime)), .groups = "drop"))

p <- ggplot(rt, aes(x = tool, y = runtime, fill = tool, colour = tool)) +
  geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
  scale_y_log10(labels = label_number(suffix = " s", big.mark = ",")) +
  scale_fill_manual(values = TOOL_COLOURS) +
  scale_colour_manual(values = TOOL_COLOURS) +
  labs(title = "Per-tool runtime — LLM family (Phase 4, real validation)",
       subtitle = "Y-axis log₁₀ scale. Each point = one of the nine real-data scenarios.",
       x = NULL, y = "Runtime (s)",
       caption = "Single-call GPTCelltype vs. multi-model consensus (CASSIA, mLLMCelltype).\nAbsolute API cost is proportional to runtime for each configuration.") +
  theme_bench()

out <- file.path(OUT_DIR, "plot_S7_llm_runtime_boxplot.png")
ggsave(out, p, width = 180, height = 140, units = "mm", dpi = 600)
cat("\nSaved:", out, "\n")
