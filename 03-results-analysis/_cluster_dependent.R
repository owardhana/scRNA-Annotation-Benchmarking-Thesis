# ============================================================
# Cluster-dependent tool decoration
# ------------------------------------------------------------
# Sourced by each phase's Plots_Analysis*.R script.
#
# Marks tools whose annotation step structurally requires a
# cluster-identity vector. In this benchmark these tools were
# supplied the ground-truth label vector as the cluster argument
# (Phases 1-3) or the per-cluster marker oracle derived on the
# full dataset (Phase 4); see Methods, "Input standardisation,
# oracle assumption, and benchmarking harness".
#
# Visual convention: dagger (U+2020) suffix on tool labels in
# axes, legends, and in-plot text annotations. Dagger is used
# in preference to asterisk to avoid collision with the
# significance-star convention (* / ** / ***) already in use
# on eta-squared bars and Nemenyi diagrams.
#
# Usage from a phase script (relative paths work from each
# results/<Phase>/ folder):
#
#   source("../_cluster_dependent.R")
#
#   # On an axis:
#   ggplot(df, aes(x = tool, ...)) + ... +
#     scale_x_discrete(labels = mark_cluster_dep) +
#     labs(..., caption = CLUSTER_DEP_CAPTION) +
#     theme_bench() + CLUSTER_DEP_THEME
#
#   # On a text label:
#   geom_text_repel(aes(label = mark_cluster_dep(tool)))
#
# The underlying `tool` column is NEVER overwritten, so joins
# and filters that key on the canonical tool name continue to
# work.
# ============================================================

CLUSTER_DEPENDENT_TOOLS <- c(
  # ----------------------------------------------------------------
  # Cell-level configurations — Supplementary Table S1
  # "Cluster-dependent" column = Yes
  # ----------------------------------------------------------------

  # Marker-based tools that consume a cluster vector
  "scType",
  "scCATCH",
  "SCSA",                # database configuration; cluster-dependent per S1
  "clustifyr_hyper",
  "clustifyr_jaccard",
  "scInfeR",

  # Correlation-based tools that consume a cluster vector
  "CIPR",
  "scmap_cluster",       # NOTE: scmap_cell is NOT cluster-dependent and is excluded

  # ----------------------------------------------------------------
  # Phase 4 database arm (P7a) — "_database" suffix variants
  # IMPORTANT: the database family is structurally heterogeneous.
  # Only the tools whose *algorithm* requires a cluster vector are
  # flagged. SCINA, scSorter, and ScInfeR are cell-by-cell marker
  # methods that, in their Phase 4 configuration, consume markers
  # sourced from a database rather than the training partition —
  # they do NOT require clusters and are therefore NOT flagged.
  # ----------------------------------------------------------------
  "scType_database",
  "scCATCH_database",
  "SCSA_database",
  "clustifyr_hyper_database",
  "clustifyr_jaccard_database",
  "scInfeR_database",
  # Intentionally excluded (cell-by-cell database-marker tools):
  #   SCINA_database, scSorter_database

  # ----------------------------------------------------------------
  # Phase 4 LLM pipelines (cluster-level input by design)
  # ----------------------------------------------------------------
  "GPTCelltype",
  "CASSIA",
  "mLLMCelltype"
)

#' Append a dagger to cluster-dependent tool names.
#'
#' Membership is an exact match against CLUSTER_DEPENDENT_TOOLS.
#' The list enumerates the cluster-dependent tools across all
#' phases of the benchmark, including the structurally
#' cluster-dependent subset of the Phase 4 database arm. Note that
#' the database family is intentionally split: scType_database /
#' scCATCH_database / SCSA_database / clustifyr_*_database output
#' labels per cluster and are flagged, while SCINA_database /
#' scSorter_database / ScInfeR_database operate cell-by-cell on
#' database markers and are NOT flagged.
#'
#' Vectorised over character / factor input. Non-matching names
#' are returned unchanged. Factor input is coerced to character;
#' callers that need a factor back should re-factor after the call.
#'
#' @param x  character or factor vector of tool names
#' @return   character vector with " †" appended to matching names
mark_cluster_dep <- function(x) {
  x_chr <- as.character(x)
  ifelse(x_chr %in% CLUSTER_DEPENDENT_TOOLS,
         paste0(x_chr, " †"),
         x_chr)
}

#' Standard plot-caption fragment explaining the dagger.
#'
#' Use as the `caption =` argument to `labs()` on any plot
#' whose tool-axis or tool-label aesthetic contains at least
#' one cluster-dependent tool.
CLUSTER_DEP_CAPTION <- paste0(
  "† Cluster-dependent tool: supplied with ground-truth labels as the cluster\n",
  "argument (Phase 1–3) or as the upstream input from which per-cluster markers\n",
  "were derived (Phase 4). See Methods, \"Input standardisation, oracle\n",
  "assumption, and benchmarking harness\"."
)

#' Theme increment for plots that use CLUSTER_DEP_CAPTION.
#'
#' The dagger caption spans 4 lines and is right-aligned by default,
#' which clips off the right edge of saved PNGs at typical widths.
#' This increment makes the caption left-aligned, lightly italicised,
#' and gives the plot enough bottom padding to render every line.
#' Add with `+` after `theme_bench()` (or your local theme).
CLUSTER_DEP_THEME <- ggplot2::theme(
  plot.caption = ggplot2::element_text(
    hjust  = 0,
    size   = 8,
    colour = "grey30",
    margin = ggplot2::margin(t = 8)
  ),
  plot.margin = ggplot2::margin(10, 12, 16, 12)
)
