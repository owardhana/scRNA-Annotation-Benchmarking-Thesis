# normalisation_pipeline.R
#
# LLM output normalisation pipeline for cell type annotation benchmarking.
#
# Maps LLM free-text cluster predictions → one of K ground-truth labels (or "unassigned")
# via a 5-stage cascade: exact match → synonym lookup → Jaro-Winkler fuzzy string →
# Cell Ontology Wu-Palmer similarity → unmappable.
#
# Produces both STRICT (WP == 1.0) and LENIENT (WP >= 0.8) mappings so that
# two kappa values can be computed per tool × dataset. The gap between them
# is itself a reported finding.
#
# Pre-registered decisions (DO NOT CHANGE):
#   - Wu-Palmer lenient threshold : 0.8
#   - Jaro-Winkler threshold      : 0.92 (applied to SIMILARITY, not distance)
#   - Unmappable outputs          → "unassigned" (not dropped, not imputed)
#   - Temperature for LLM calls   : 0  (set in calling tool files)
#
# Source this file from the benchmarking root:
#   source("LLM-based/normalisation_pipeline.R", local = TRUE)

# ---------------------------------------------------------------------------
# Package availability check (warn, do not error)
# ---------------------------------------------------------------------------
.required_pkgs <- c(
  "ontologyIndex",      # get_ancestors(), get_OBO()
  "ontologySimilarity", # loaded for completeness — Wu-Palmer is implemented manually
  "ontoProc",           # getOnto("cellOnto")
  "stringdist",         # Jaro-Winkler via stringdist(..., method = "jw")
  "stringr",            # str_squish(), str_replace_all()
  "dplyr",              # data manipulation
  "irr",                # kappa2() for Cohen's kappa
  "caret",              # confusionMatrix() for macro F1
  "MLmetrics"           # MCC fallback
)
.missing_pkgs <- .required_pkgs[!sapply(.required_pkgs, requireNamespace, quietly = TRUE)]
if (length(.missing_pkgs) > 0) {
  warning(
    "normalisation_pipeline.R: the following packages are not installed and must be ",
    "installed before using all pipeline functions:\n  ",
    paste(.missing_pkgs, collapse = ", ")
  )
}
rm(.required_pkgs, .missing_pkgs)


# ===========================================================================
# FUNCTION 1: load_cell_ontology
# ===========================================================================

#' Load the Cell Ontology (CL) as an ontology_index object
#'
#' @return An \code{ontology_index} object with names/ancestors for CL terms.
#'
#' @details
#' Tries \code{ontoProc::getOnto("cellOnto")} first. If that fails, falls back
#' to reading \code{resources/cl.obo} via \code{ontologyIndex::get_OBO()}.
#' Prints which method succeeded.
load_cell_ontology <- function() {
  cl_onto <- tryCatch({
    if (!requireNamespace("ontoProc", quietly = TRUE)) stop("ontoProc not available")
    onto <- ontoProc::getOnto("cellOnto")
    cat("[normalisation_pipeline] Cell Ontology loaded via ontoProc::getOnto()\n")
    onto
  }, error = function(e) {
    cat("[normalisation_pipeline] ontoProc failed (", conditionMessage(e),
        ") — falling back to get_OBO()\n", sep = "")
    if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
      stop("ontologyIndex package is required as fallback for load_cell_ontology()")
    }
    obo_path <- "resources/cl.obo"
    if (!file.exists(obo_path)) {
      stop(
        "Fallback OBO file not found at '", obo_path, "'. ",
        "Either install ontoProc or place cl.obo in resources/"
      )
    }
    onto <- ontologyIndex::get_OBO(
      obo_path,
      propagate_relationships = c("is_a", "part_of")
    )
    cat("[normalisation_pipeline] Cell Ontology loaded via ontologyIndex::get_OBO()\n")
    onto
  })
  cl_onto
}


# ===========================================================================
# FUNCTION 2: normalise_text
# ===========================================================================

#' Normalise a raw LLM output string for downstream matching
#'
#' @param raw_text A single character string (LLM output for one cluster).
#'
#' @return Cleaned lowercase string, or \code{NA} if input is NA or empty.
#'
#' @details
#' Steps applied in order:
#' \enumerate{
#'   \item Return NA for NA or empty input.
#'   \item Lowercase and trim leading/trailing whitespace.
#'   \item Strip common LLM preambles.
#'   \item Remove trailing period.
#'   \item Replace "/" with " or " (e.g. "CD4/CD8" → "CD4 or CD8").
#'   \item Collapse multiple spaces to single space.
#' }
normalise_text <- function(raw_text) {
  if (is.na(raw_text) || !nzchar(trimws(raw_text))) return(NA_character_)

  x <- tolower(trimws(raw_text))

  # Remove common LLM preambles (order matters: longer patterns first)
  preambles <- c(
    "based on the markers, ",
    "i would classify this as ",
    "this is a ",
    "cell type: ",
    "the "
  )
  for (p in preambles) {
    if (startsWith(x, p)) {
      x <- substr(x, nchar(p) + 1L, nchar(x))
      break   # only strip one preamble
    }
  }

  # Remove trailing period
  x <- sub("\\.$", "", x)

  # Replace "/" with " or "
  x <- gsub("/", " or ", x, fixed = TRUE)

  # Collapse multiple spaces
  x <- trimws(gsub("\\s+", " ", x))

  if (!nzchar(x)) NA_character_ else x
}


# ===========================================================================
# FUNCTION 3: build_synonym_table
# ===========================================================================

#' Build a synonym look-up table filtered to the ground-truth labels of a dataset
#'
#' @param ground_truth_labels Character vector of all unique ground-truth cell
#'   type labels present in this dataset.
#'
#' @return A named list mapping variant strings → canonical strings (both
#'   lowercase). Only entries whose canonical form appears in
#'   \code{tolower(ground_truth_labels)} are returned.
#'
#' @details
#' The hard-coded synonym map covers at minimum 24 canonical cell-type groups.
#' Variants are common abbreviations and wordings found in LLM outputs.
build_synonym_table <- function(ground_truth_labels) {
  # Hard-coded synonym map: variant → canonical (all lowercase)
  raw_map <- list(
    # T cell
    "t cell"                     = "t cell",
    "t-cell"                     = "t cell",
    "t lymphocyte"               = "t cell",
    "t-lymphocyte"               = "t cell",
    "t cells"                    = "t cell",

    # B cell
    "b cell"                     = "b cell",
    "b-cell"                     = "b cell",
    "b lymphocyte"               = "b cell",
    "b cells"                    = "b cell",

    # NK cell
    "nk cell"                    = "nk cell",
    "nk cells"                   = "nk cell",
    "natural killer cell"        = "nk cell",
    "natural killer cells"       = "nk cell",
    "natural killer"             = "nk cell",

    # Monocyte (generic)
    "monocyte"                   = "monocyte",
    "monocytes"                  = "monocyte",

    # Classical monocyte
    "classical monocyte"         = "classical monocyte",
    "cd14+ monocyte"             = "classical monocyte",
    "cd14 monocyte"              = "classical monocyte",
    "cd14+ monocytes"            = "classical monocyte",
    "cd14 monocytes"             = "classical monocyte",
    "classical monocytes"        = "classical monocyte",

    # Non-classical monocyte
    "non-classical monocyte"     = "non-classical monocyte",
    "nonclassical monocyte"      = "non-classical monocyte",
    "cd16+ monocyte"             = "non-classical monocyte",
    "cd16 monocyte"              = "non-classical monocyte",
    "cd16+ monocytes"            = "non-classical monocyte",
    "non-classical monocytes"    = "non-classical monocyte",

    # Dendritic cell
    "dendritic cell"             = "dendritic cell",
    "dendritic cells"            = "dendritic cell",
    "dc"                         = "dendritic cell",
    "myeloid dendritic cell"     = "dendritic cell",

    # Plasmacytoid dendritic cell
    "plasmacytoid dendritic cell"  = "plasmacytoid dendritic cell",
    "plasmacytoid dendritic cells" = "plasmacytoid dendritic cell",
    "pdc"                          = "plasmacytoid dendritic cell",
    "pdcs"                         = "plasmacytoid dendritic cell",
    "plasmacytoid dc"              = "plasmacytoid dendritic cell",

    # Macrophage
    "macrophage"                 = "macrophage",
    "macrophages"                = "macrophage",

    # Neutrophil
    "neutrophil"                 = "neutrophil",
    "neutrophils"                = "neutrophil",

    # CD8 T cell
    "cd8 t cell"                 = "cd8 t cell",
    "cd8+ t cell"                = "cd8 t cell",
    "cd8+ t cells"               = "cd8 t cell",
    "cd8 t cells"                = "cd8 t cell",
    "cytotoxic t cell"           = "cd8 t cell",
    "cytotoxic t cells"          = "cd8 t cell",
    "cytotoxic t lymphocyte"     = "cd8 t cell",
    "ctl"                        = "cd8 t cell",

    # CD4 T cell
    "cd4 t cell"                 = "cd4 t cell",
    "cd4+ t cell"                = "cd4 t cell",
    "cd4+ t cells"               = "cd4 t cell",
    "cd4 t cells"                = "cd4 t cell",
    "helper t cell"              = "cd4 t cell",
    "helper t cells"             = "cd4 t cell",
    "th cell"                    = "cd4 t cell",

    # Regulatory T cell
    "regulatory t cell"          = "regulatory t cell",
    "regulatory t cells"         = "regulatory t cell",
    "treg"                       = "regulatory t cell",
    "tregs"                      = "regulatory t cell",
    "t regulatory cell"          = "regulatory t cell",
    "foxp3+ t cell"              = "regulatory t cell",

    # Naive CD4 T cell
    "naive cd4 t cell"           = "naive cd4 t cell",
    "naive cd4+ t cell"          = "naive cd4 t cell",
    "naive cd4 t cells"          = "naive cd4 t cell",

    # Memory CD4 T cell
    "memory cd4 t cell"          = "memory cd4 t cell",
    "memory cd4+ t cell"         = "memory cd4 t cell",
    "memory cd4 t cells"         = "memory cd4 t cell",

    # Naive CD8 T cell
    "naive cd8 t cell"           = "naive cd8 t cell",
    "naive cd8+ t cell"          = "naive cd8 t cell",
    "naive cd8 t cells"          = "naive cd8 t cell",

    # Memory CD8 T cell
    "memory cd8 t cell"          = "memory cd8 t cell",
    "memory cd8+ t cell"         = "memory cd8 t cell",
    "memory cd8 t cells"         = "memory cd8 t cell",

    # Natural killer T cell
    "natural killer t cell"      = "natural killer t cell",
    "nkt cell"                   = "natural killer t cell",
    "nkt cells"                  = "natural killer t cell",
    "nk t cell"                  = "natural killer t cell",

    # Plasma cell
    "plasma cell"                = "plasma cell",
    "plasma cells"               = "plasma cell",
    "plasmablast"                = "plasma cell",
    "antibody-secreting cell"    = "plasma cell",

    # Platelet / thrombocyte
    "platelet"                   = "platelet",
    "platelets"                  = "platelet",
    "thrombocyte"                = "platelet",

    # Erythrocyte / red blood cell
    "erythrocyte"                = "erythrocyte",
    "erythrocytes"               = "erythrocyte",
    "red blood cell"             = "erythrocyte",
    "rbc"                        = "erythrocyte",

    # Mast cell
    "mast cell"                  = "mast cell",
    "mast cells"                 = "mast cell",

    # Basophil
    "basophil"                   = "basophil",
    "basophils"                  = "basophil",

    # Eosinophil
    "eosinophil"                 = "eosinophil",
    "eosinophils"                = "eosinophil"
  )

  gt_lower <- tolower(ground_truth_labels)
  # Keep only entries whose canonical appears in the dataset's GT labels
  filtered <- Filter(function(canonical) canonical %in% gt_lower, raw_map)
  filtered
}


# ===========================================================================
# FUNCTION 4: text_to_cl_term
# ===========================================================================

#' Map a normalised text label to a Cell Ontology term ID
#'
#' @param text_label Normalised lowercase string (output of \code{normalise_text()}).
#' @param cl_onto An \code{ontology_index} object (output of \code{load_cell_ontology()}).
#'
#' @return A CL term ID string (e.g. \code{"CL:0000084"}) or \code{NA} if no
#'   match is found. Returns \code{NA} gracefully if \code{text_label} is \code{NA}.
#'
#' @details
#' Step 1: exact case-insensitive match against \code{cl_onto$name}.
#' Step 2: partial match — \code{text_label} appears as a substring in a CL label.
#' Returns the first hit in both steps.
text_to_cl_term <- function(text_label, cl_onto) {
  if (is.na(text_label) || !nzchar(text_label)) return(NA_character_)

  cl_names_lower <- tolower(cl_onto$name)

  # Step 1: exact match
  exact_idx <- which(cl_names_lower == text_label)
  if (length(exact_idx) > 0L) return(names(cl_onto$name)[exact_idx[1L]])
  
  # Step 2: whole-word match
  word_idx <- which(grepl(paste0("\\b", text_label, "\\b"), 
                          cl_names_lower, perl = TRUE))
  if (length(word_idx) > 0L) {
    label_lengths <- nchar(cl_names_lower[word_idx])
    best <- word_idx[which.min(label_lengths)]
    return(names(cl_onto$name)[best])
  }

  # Step 3: partial match — rank by label length ascending to prefer the
  # shortest (most specific) CL term that contains text_label as a substring
  partial_idx <- which(grepl(text_label, cl_names_lower, fixed = TRUE))
  if (length(partial_idx) > 0L) {
    label_lengths <- nchar(cl_names_lower[partial_idx])
    best_partial  <- partial_idx[which.min(label_lengths)]
    return(names(cl_onto$name)[best_partial])
  }

  NA_character_
}


# ===========================================================================
# FUNCTION 5: wu_palmer_similarity
# ===========================================================================

#' Compute Wu-Palmer semantic similarity between two Cell Ontology terms
#'
#' @param term_a A CL term ID string (e.g. \code{"CL:0000084"}).
#' @param term_b A CL term ID string.
#' @param cl_onto An \code{ontology_index} object.
#'
#' @return Numeric in [0, 1]. Returns 1.0 if \code{term_a == term_b}, 0.0 if
#'   either is \code{NA} or if the terms share no common ancestor.
#'
#' @details
#' Wu-Palmer similarity is implemented manually as:
#' \deqn{WP(A, B) = \frac{2 \cdot \text{depth}(LCS)}{\text{depth}(A) + \text{depth}(B)}}
#' where \eqn{\text{depth}(X) = |\text{ancestors}(X)|}  and
#' LCS is the Least Common Subsumer (common ancestor with maximum depth).
#'
#' Uses \code{ontologyIndex::get_ancestors()} — NOT deprecated \code{ancestors()}.
#'Note: get_ancestors() includes the term itself, so depth(X) = 
# number of ancestors + 1. This is consistent across both terms
# and does not affect the relative ordering of WP scores.
wu_palmer_similarity <- function(term_a, term_b, cl_onto) {
  if (is.na(term_a) || is.na(term_b)) return(0.0)
  if (term_a == term_b) return(1.0)

  anc_a <- ontologyIndex::get_ancestors(cl_onto, term_a)
  anc_b <- ontologyIndex::get_ancestors(cl_onto, term_b)

  depth_a <- length(anc_a)
  depth_b <- length(anc_b)

  common <- intersect(anc_a, anc_b)
  if (length(common) == 0L) return(0.0)

  # Depth of each common ancestor = number of its own ancestors
  common_depths <- sapply(common, function(x) length(ontologyIndex::get_ancestors(cl_onto, x)))
  depth_lcs <- max(common_depths)

  denom <- depth_a + depth_b
  if (denom == 0L) return(0.0)

  (2.0 * depth_lcs) / denom
}


# ===========================================================================
# FUNCTION 6: map_single_prediction
# ===========================================================================

#' Map one normalised cluster prediction to a ground-truth label
#'
#' @param normalised_pred Normalised lowercase string from \code{normalise_text()}.
#' @param synonym_lookup Named list from \code{build_synonym_table()}.
#' @param ground_truth_labels Character vector of all unique GT labels in this dataset.
#' @param cl_onto An \code{ontology_index} object.
#' @param wp_threshold Numeric; Wu-Palmer threshold for lenient match (default 0.8).
#' @param jw_threshold Numeric; Jaro-Winkler SIMILARITY threshold (default 0.92).
#'
#' @return Named list with fields:
#' \describe{
#'   \item{strict}{Matched GT label or \code{"unassigned"}.}
#'   \item{lenient}{Matched GT label or \code{"unassigned"}.}
#'   \item{wp_sim}{Best Wu-Palmer similarity achieved, or \code{NA} if not reached.}
#'   \item{mapping_method}{One of \code{"na_input"}, \code{"exact"}, \code{"synonym"},
#'     \code{"fuzzy_string"}, \code{"ontology_wp"}, \code{"unmappable"}.}
#' }
map_single_prediction <- function(normalised_pred,
                                  synonym_lookup,
                                  ground_truth_labels,
                                  cl_onto,
                                  wp_threshold = 0.8,
                                  jw_threshold = 0.92) {
  unassigned_result <- function(method, wp = NA_real_) {
    list(strict = "Unknown", lenient = "Unknown",
         wp_sim = wp, mapping_method = method)
  }
  match_result <- function(gt_label, method, wp = 1.0) {
    list(strict = gt_label, lenient = gt_label,
         wp_sim = wp, mapping_method = method)
  }

  # Stage 0: NA / empty input
  if (is.na(normalised_pred) || !nzchar(normalised_pred)) {
    return(unassigned_result("na_input"))
  }

  gt_lower <- tolower(ground_truth_labels)

  # Stage 1: Exact match
  exact_idx <- match(normalised_pred, gt_lower)
  if (!is.na(exact_idx)) {
    return(match_result(ground_truth_labels[exact_idx], "exact"))
  }

  # Stage 2: Synonym lookup
  canonical <- synonym_lookup[[normalised_pred]]
  if (!is.null(canonical) && canonical %in% gt_lower) {
    gt_match <- ground_truth_labels[match(canonical, gt_lower)]
    return(match_result(gt_match, "synonym"))
  }

  # Stage 3: Jaro-Winkler fuzzy string
  # stringdist() returns DISTANCE; convert to similarity = 1 - distance
  jw_distances <- stringdist::stringdist(normalised_pred, gt_lower, method = "jw")
  jw_sims <- 1.0 - jw_distances
  best_jw_idx <- which.max(jw_sims)
  if (length(best_jw_idx) > 0L && jw_sims[best_jw_idx] >= jw_threshold) {
    return(match_result(ground_truth_labels[best_jw_idx], "fuzzy_string"))
  }

  # Stage 4: Cell Ontology Wu-Palmer
  if (is.null(cl_onto)) return(unassigned_result("unmappable"))
  pred_cl <- text_to_cl_term(normalised_pred, cl_onto)
  if (!is.na(pred_cl)) {
    gt_cl_terms <- sapply(gt_lower, text_to_cl_term, cl_onto = cl_onto,
                          USE.NAMES = FALSE)
    valid_idx <- which(!is.na(gt_cl_terms))
    if (length(valid_idx) > 0L) {
      wp_sims <- sapply(valid_idx, function(i) {
        wu_palmer_similarity(pred_cl, gt_cl_terms[i], cl_onto)
      })
      best_i  <- which.max(wp_sims)
      best_wp <- wp_sims[best_i]
      best_gt <- ground_truth_labels[valid_idx[best_i]]
      
      # strict requires exact CL node match (wp == 1.0).
      # In practice this is unreachable from Stage 4 since exact/synonym
      # matching would have resolved it earlier. Stage 4 strict labels
      # will always be "unassigned", which is the intended behaviour —
      # it penalises granularity mismatches in the strict kappa metric.
      strict_label  <- if (best_wp == 1.0) best_gt else "Unknown"
      lenient_label <- if (best_wp >= wp_threshold) best_gt else "Unknown"

      return(list(
        strict         = strict_label,
        lenient        = lenient_label,
        wp_sim         = best_wp,
        mapping_method = "ontology_wp"
      ))
    }
  }

  # Stage 5: Unmappable
  unassigned_result("unmappable")
}


# ===========================================================================
# FUNCTION 7: run_normalisation_pipeline
# ===========================================================================

#' Run the full normalisation pipeline over all cluster predictions
#'
#' @param llm_raw_df A data.frame with columns \code{cluster} and \code{raw_pred}.
#' @param ground_truth_labels Character vector of all unique GT labels in this dataset.
#' @param cl_onto An \code{ontology_index} object from \code{load_cell_ontology()}.
#' @param wp_threshold Wu-Palmer lenient threshold (default 0.8).
#' @param jw_threshold Jaro-Winkler similarity threshold (default 0.92).
#'
#' @return A data.frame with columns:
#'   \code{cluster}, \code{raw_pred}, \code{normalised_pred},
#'   \code{strict}, \code{lenient}, \code{wp_sim}, \code{mapping_method}.
#'
#' @details
#' Builds the synonym lookup once, then applies \code{normalise_text()} and
#' \code{map_single_prediction()} to each row. Prints a summary of how many
#' clusters fell into each \code{mapping_method} category.
run_normalisation_pipeline <- function(llm_raw_df,
                                       ground_truth_labels,
                                       cl_onto,
                                       wp_threshold = 0.8,
                                       jw_threshold = 0.92) {
  stopifnot(all(c("cluster", "raw_pred") %in% names(llm_raw_df)))

  synonym_lookup <- build_synonym_table(ground_truth_labels)

  results <- lapply(seq_len(nrow(llm_raw_df)), function(i) {
    norm_pred <- normalise_text(llm_raw_df$raw_pred[i])
    mapped    <- map_single_prediction(
      norm_pred, synonym_lookup, ground_truth_labels, cl_onto,
      wp_threshold, jw_threshold
    )
    list(
      cluster        = llm_raw_df$cluster[i],
      raw_pred       = llm_raw_df$raw_pred[i],
      normalised_pred = norm_pred,
      strict         = mapped$strict,
      lenient        = mapped$lenient,
      wp_sim         = mapped$wp_sim,
      mapping_method = mapped$mapping_method
    )
  })

  out_df <- data.frame(
    cluster         = sapply(results, `[[`, "cluster"),
    raw_pred        = sapply(results, `[[`, "raw_pred"),
    normalised_pred = sapply(results, `[[`, "normalised_pred"),
    strict          = sapply(results, `[[`, "strict"),
    lenient         = sapply(results, `[[`, "lenient"),
    wp_sim          = as.numeric(sapply(results, `[[`, "wp_sim")),
    mapping_method  = sapply(results, `[[`, "mapping_method"),
    stringsAsFactors = FALSE
  )

  # Summary
  method_counts <- table(out_df$mapping_method)
  cat("\n[normalisation_pipeline] Mapping method summary (", nrow(out_df), " clusters):\n", sep = "")
  print(method_counts)
  cat("\n")

  out_df
}


# ===========================================================================
# FUNCTION 8: build_confusion_and_compute_metrics
# ===========================================================================

#' Build confusion matrix and compute strict/lenient metrics from normalised predictions
#'
#' @param test_cell_df A data.frame with columns \code{cell_id}, \code{cluster},
#'   \code{true_label} — one row per test cell.
#' @param normalised_df Output of \code{run_normalisation_pipeline()}.
#' @param ground_truth_labels Character vector of all unique GT labels in this dataset.
#'
#' @return A named list:
#' \describe{
#'   \item{strict}{list(kappa, macro_f1, accuracy, mcc, unassigned_rate)}
#'   \item{lenient}{list(kappa, macro_f1, accuracy, mcc, unassigned_rate)}
#'   \item{hallucination_rate}{Proportion of CLUSTERS with method == "unmappable".}
#'   \item{kappa_gap}{lenient$kappa - strict$kappa.}
#'   \item{per_cell_df}{Full joined data.frame for inspection.}
#' }
build_confusion_and_compute_metrics <- function(test_cell_df,
                                                normalised_df,
                                                ground_truth_labels) {
  stopifnot(all(c("cell_id", "cluster", "true_label") %in% names(test_cell_df)))
  stopifnot(all(c("cluster", "strict", "lenient", "mapping_method") %in% names(normalised_df)))

  # Join on cluster
  joined <- merge(test_cell_df, normalised_df[, c("cluster", "strict", "lenient", "mapping_method")],
                  by = "cluster", all.x = TRUE)
  joined$strict[is.na(joined$strict)]   <- "Unknown"
  joined$lenient[is.na(joined$lenient)] <- "Unknown"

  # All possible levels: GT labels + "Unknown"
  all_levels <- c(sort(ground_truth_labels), "Unknown")

  .compute_set <- function(pred_vec, true_vec) {
    pred_f <- factor(pred_vec, levels = all_levels)
    true_f <- factor(true_vec, levels = all_levels)

    # Cohen's kappa
    kappa_val <- tryCatch({
      kappa_df <- data.frame(r1 = true_f, r2 = pred_f)
      irr::kappa2(kappa_df)$value
    }, error = function(e) NA_real_)

    # confusionMatrix (suppress warnings from missing levels)
    cm <- tryCatch(
      suppressWarnings(caret::confusionMatrix(pred_f, true_f)),
      error = function(e) NULL
    )

    accuracy <- if (!is.null(cm)) unname(cm$overall["Accuracy"]) else NA_real_

    macro_f1 <- if (!is.null(cm)) {
      f1_per_class <- cm$byClass[, "F1"]
      mean(f1_per_class, na.rm = TRUE)
    } else NA_real_

    # MCC — Gorodkin (2004) multiclass formula
    mcc_val <- tryCatch({
      conf_mat <- table(true = true_f, pred = pred_f)
      c_val <- sum(diag(conf_mat))
      s_val <- sum(conf_mat)
      pk <- colSums(conf_mat)   # predicted counts per class
      tk <- rowSums(conf_mat)   # true counts per class
      numer <- c_val * s_val - sum(pk * tk)
      denom <- sqrt((s_val^2 - sum(pk^2)) * (s_val^2 - sum(tk^2)))
      if (is.na(denom) || denom == 0) 0.0 else numer / denom
    }, error = function(e) NA_real_)

    unassigned_rate <- mean(pred_vec == "Unknown", na.rm = TRUE)

    list(
      kappa          = kappa_val,
      macro_f1       = macro_f1,
      accuracy       = accuracy,
      mcc            = mcc_val,
      unassigned_rate = unassigned_rate
    )
  }

  strict_metrics  <- .compute_set(joined$strict,  joined$true_label)
  lenient_metrics <- .compute_set(joined$lenient, joined$true_label)

  # Hallucination rate: proportion of CLUSTERS (not cells) that are unmappable
  hallucination_rate <- mean(normalised_df$mapping_method == "unmappable", na.rm = TRUE)

  kappa_gap <- tryCatch(lenient_metrics$kappa - strict_metrics$kappa, error = function(e) NA_real_)

  list(
    strict           = strict_metrics,
    lenient          = lenient_metrics,
    hallucination_rate = hallucination_rate,
    kappa_gap        = kappa_gap,
    per_cell_df      = joined
  )
}