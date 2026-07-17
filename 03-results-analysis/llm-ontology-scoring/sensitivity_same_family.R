# ============================================================
# Phase 4 LLM-panel scoring — same-family sensitivity analysis
# ============================================================
# Companion to score_cells.py. Reads every *_cluster_map.tsv in
# ./with_score/, computes the full-panel mean (score_1, score_2,
# score_3 averaged) and the leave-one-family-out sensitivity means
# per tool, then writes results/Sensitivity_SameFamily.csv.
#
# Scorer -> family mapping (per score_cells.py):
#   score_1 = google/gemini-3-flash-preview   -> family: Gemini
#   score_2 = anthropic/claude-haiku-4-5      -> family: Claude
#   score_3 = openai/gpt-5.4-nano             -> family: GPT
#
# Tool -> family mapping (which model families each tool emits from):
#   GPTCelltype  -> {GPT}
#   CASSIA       -> {GPT, Claude}
#   mLLMCelltype -> {GPT, Claude, Gemini}
#
# Sensitivity scheme:
#   * Tools with at least one cross-family scorer remaining
#     (GPTCelltype, CASSIA): strict same-family removal.
#   * mLLMCelltype: all three scorers are same-family -> strict
#     removal is undefined. Reported under leave-one-family-out:
#     three means, each dropping one scorer singly.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

SCORE_DIR <- "./with_score"
OUT_DIR   <- "./results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Scorer -> family
SCORER_FAMILY <- c(score_1 = "Gemini", score_2 = "Claude", score_3 = "GPT")

# Tool -> set of families it emits from
TOOL_FAMILIES <- list(
  GPTCelltype  = c("GPT"),
  CASSIA       = c("GPT", "Claude"),
  mLLMCelltype = c("GPT", "Claude", "Gemini")
)

# ---- 1. Discover and load scored TSVs ----------------------

discover_scored_files <- function(dir = SCORE_DIR) {
  files <- list.files(dir, pattern = "_cluster_map\\.tsv$", full.names = TRUE)
  if (length(files) == 0L)
    stop("No scored TSVs found in ", dir,
         ". Run score_cells.py first.", call. = FALSE)
  tibble::tibble(
    path    = files,
    fname   = basename(files),
    dataset = str_replace(fname, "_(GPTCelltype|CASSIA|mLLMCelltype)_cluster_map\\.tsv$", ""),
    tool    = str_extract(fname, "GPTCelltype|CASSIA|mLLMCelltype")
  )
}

load_scored <- function(row) {
  read_tsv(row$path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(dataset = row$dataset, tool = row$tool)
}

# ---- 2. Score aggregation ----------------------------------

#' Mean over a subset of scorer columns, ignoring NA per row.
#' Returns NA if all selected scorers are NA for a row.
row_panel_mean <- function(df, cols) {
  m <- as.matrix(df[, cols, drop = FALSE])
  rowMeans(m, na.rm = TRUE)
}

#' Per-tool aggregate: arithmetic mean of per-row panel means,
#' across all clusters and datasets for that tool.
agg_tool_panel <- function(df, cols) {
  per_row <- row_panel_mean(df, cols)
  per_row <- per_row[is.finite(per_row)]
  if (length(per_row) == 0L) NA_real_
  else mean(per_row)
}

# ---- 3. Same-family identification per tool ----------------

same_family_scorers <- function(tool) {
  fams <- TOOL_FAMILIES[[tool]]
  names(SCORER_FAMILY)[SCORER_FAMILY %in% fams]
}

cross_family_scorers <- function(tool) {
  fams <- TOOL_FAMILIES[[tool]]
  names(SCORER_FAMILY)[!SCORER_FAMILY %in% fams]
}

# ---- 4. Per-tool sensitivity row builder -------------------

per_tool_sensitivity <- function(tool_df, tool) {
  full_cols  <- c("score_1", "score_2", "score_3")
  full_mean  <- agg_tool_panel(tool_df, full_cols)

  same_fam <- same_family_scorers(tool)
  cross_fam <- cross_family_scorers(tool)

  # Strict same-family removal (only defined when >=1 cross-family scorer remains)
  if (length(cross_fam) >= 1L) {
    strict_mean <- agg_tool_panel(tool_df, cross_fam)
    strict_label <- paste(SCORER_FAMILY[same_fam], collapse = " + ")
    n_strict <- length(cross_fam)
  } else {
    strict_mean  <- NA_real_
    strict_label <- "undefined (all scorers same-family)"
    n_strict     <- 0L
  }

  # Leave-one-family-out: one row per scorer dropped
  loo_rows <- purrr::map_dfr(full_cols, function(drop_col) {
    keep_cols <- setdiff(full_cols, drop_col)
    tibble::tibble(
      tool             = tool,
      scheme           = "leave-one-family-out",
      dropped_scorer   = drop_col,
      dropped_family   = SCORER_FAMILY[[drop_col]],
      n_scorers_used   = length(keep_cols),
      mean_score       = agg_tool_panel(tool_df, keep_cols)
    )
  })

  summary_rows <- tibble::tibble(
    tool             = tool,
    scheme           = c("full_panel", "strict_same_family_removed"),
    dropped_scorer   = c(NA, paste(same_fam, collapse = ",")),
    dropped_family   = c(NA, strict_label),
    n_scorers_used   = c(3L, n_strict),
    mean_score       = c(full_mean, strict_mean)
  )

  dplyr::bind_rows(summary_rows, loo_rows)
}

# ---- 5. Main pipeline --------------------------------------

main <- function() {
  scored_index <- discover_scored_files()
  message("Found ", nrow(scored_index), " scored TSV(s) across ",
          dplyr::n_distinct(scored_index$tool), " tool(s) and ",
          dplyr::n_distinct(scored_index$dataset), " dataset(s).")

  all_scored <- scored_index %>%
    purrr::pmap_dfr(function(...) load_scored(tibble::tibble(...)))

  tools_present <- intersect(names(TOOL_FAMILIES), unique(all_scored$tool))

  sensitivity_df <- purrr::map_dfr(tools_present, function(t) {
    per_tool_sensitivity(dplyr::filter(all_scored, tool == t), t)
  }) %>%
    dplyr::mutate(mean_score = round(mean_score, 4))

  # Rank under each scheme (1 = best)
  ranks_df <- sensitivity_df %>%
    dplyr::filter(!is.na(mean_score)) %>%
    dplyr::group_by(scheme, dropped_scorer) %>%
    dplyr::mutate(rank = dplyr::dense_rank(-mean_score)) %>%
    dplyr::ungroup()

  out_path <- file.path(OUT_DIR, "Sensitivity_SameFamily.csv")
  readr::write_csv(ranks_df, out_path)
  message("Wrote: ", out_path)

  # Pretty console summary
  cat("\n=== Same-family sensitivity summary ===\n")
  ranks_df %>%
    dplyr::filter(scheme %in% c("full_panel", "strict_same_family_removed")) %>%
    dplyr::select(tool, scheme, mean_score, rank) %>%
    print(n = Inf)

  cat("\n=== Leave-one-family-out (per scorer dropped) ===\n")
  ranks_df %>%
    dplyr::filter(scheme == "leave-one-family-out") %>%
    dplyr::select(tool, dropped_family, mean_score, rank) %>%
    print(n = Inf)

  invisible(ranks_df)
}

if (sys.nframe() == 0L) main()
