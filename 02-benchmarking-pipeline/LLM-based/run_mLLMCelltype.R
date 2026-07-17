# run_mLLMCelltype.R
#
# mLLMCelltype 3-model collaboration variant:
#   - claude-sonnet-4-6  (Anthropic, via OpenRouter)
#   - gpt-5.4            (OpenAI, direct or via OpenRouter)
#   - google/gemini-3.1-pro-preview (Google, via OpenRouter)
#
# Uses interactive_consensus_annotation() for multi-model discussion,
# then maps cluster-level LLM text → GT labels via normalisation_pipeline.R,
# yielding both STRICT and LENIENT predictions / kappa values.

run_mLLMCelltype_function <- function(seurat_train, seurat_test, markers) {
  # Load required libraries with error checking
  required_packages <- c("mLLMCelltype", "Seurat")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "not available. Please install it first."))
    }
  }
  library(mLLMCelltype)
  library(Seurat)

  # Use centralized prepare_markers() helper for uniform filtering
  source("benchmarking_helpers.R", local = TRUE)
  markers <- prepare_markers(markers)

  # Load normalisation pipeline
  source("LLM-based/normalisation_pipeline.R", local = TRUE)

  # Default return function for error handling
  default_return <- function() {
    n <- ncol(seurat_test)
    return(list(
      predictions = rep("Unknown", n),
      lenient_predictions = rep("Unknown", n),
      true_labels = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores = rep(0, n),
      cell_ids = colnames(seurat_test),
      strict_kappa = NA_real_,
      lenient_kappa = NA_real_,
      kappa_gap = NA_real_,
      hallucination_rate = NA_real_,
      runtime_secs = NA_real_,
      peak_system_memory_mb = NA_real_
    ))
  }

  # ---- Configure API keys -------------------------------------------------
  api_keys <- list(
    openrouter = Sys.getenv("OPENROUTER_API_KEY"),
    openai     = Sys.getenv("OPENAI_API_KEY"),
    gemini     = Sys.getenv("GEMINI_API_KEY")
  )

  if (nchar(api_keys$openrouter) == 0 && nchar(api_keys$openai) == 0) {
    warning("mLLMCelltype: no API keys found (OPENROUTER_API_KEY or OPENAI_API_KEY required)")
    return(default_return())
  }

  # ---- 3-model ensemble ---------------------------------------------------
  models <- c(
    "anthropic/claude-sonnet-4-6",
    "gpt-5.4",
    "google/gemini-3.1-pro-preview"
  )

  # ---- Map cluster IDs ----------------------------------------------------
  # In CV mode, FindAllMarkers runs with Ground_Truth_Celltype identities, so
  # markers$cluster holds cell type names that must be remapped to cluster IDs.
  # In full-dataset mode, FindAllMarkers runs with Ground_Truth_Cluster identities,
  # so markers$cluster already holds cluster IDs — skip the remap.
  if (!"Ground_Truth_Cluster" %in% colnames(seurat_test@meta.data)) {
    warning("mLLMCelltype: Ground_Truth_Cluster column missing from seurat_test — cannot build cell data frame")
    return(default_return())
  }
  if ("Ground_Truth_Cluster" %in% colnames(seurat_train@meta.data)) {
    mapping_df <- unique(seurat_train@meta.data[, c("Ground_Truth_Celltype", "Ground_Truth_Cluster")])
    celltype_to_cluster_map <- setNames(mapping_df$Ground_Truth_Cluster, mapping_df$Ground_Truth_Celltype)
    # Only remap if markers are keyed by Ground_Truth_Celltype (CV mode)
    if (any(as.character(markers$cluster) %in% mapping_df$Ground_Truth_Celltype)) {
      markers$cluster <- celltype_to_cluster_map[as.character(markers$cluster)]
      na_marker_rows <- sum(is.na(markers$cluster))
      if (na_marker_rows > 0) {
        warning("mLLMCelltype: ", na_marker_rows, " marker rows have no cluster mapping and will be removed")
        markers <- markers[!is.na(markers$cluster), ]
      }
      if (nrow(markers) == 0) {
        warning("mLLMCelltype: no markers remain after cluster remapping")
        return(default_return())
      }
    }
  }

  # Track peak memory usage
  runtime_secs <- NA_real_
  peak_system_memory_mb <- NA_real_

  # ---- Run consensus annotation (only this is memory-tracked) -------------
  consensus_results <- NULL

  run_consensus <- function() {
    tryCatch(
      {
        cat("mLLMCelltype: calling interactive_consensus_annotation with 3 models...\n")
        mLLMCelltype::interactive_consensus_annotation(
          input = markers,
          tissue_name = "generic",
          top_gene_count = 20,
          models = models,
          api_keys = api_keys,
          controversy_threshold = 0.7,
          entropy_threshold = 1.0,
          max_discussion_rounds = 3,
          consensus_check_model = NULL,
          force_rerun = TRUE,
          use_cache = FALSE
        )
      },
      error = function(e) {
        warning("mLLMCelltype consensus annotation failed: ", e$message)
        NULL
      }
    )
  }

  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
    start_time <- Sys.time()
    consensus_results <- run_consensus()
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM({
      consensus_results <- run_consensus()
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  if (is.null(consensus_results)) {
    warning("mLLMCelltype: consensus failed — returning default predictions")
    return(default_return())
  }

  # ---- Extract cluster-level raw annotations ------------------------------
  cluster_to_celltype_map <- consensus_results$final_annotations
  unique_cluster_ids <- sort(unique(markers$cluster))

  llm_raw_df <- data.frame(
    cluster = as.character(unique_cluster_ids),
    raw_pred = as.character(cluster_to_celltype_map[as.character(unique_cluster_ids)]),
    stringsAsFactors = FALSE
  )
  # Replace NA raw_pred (clusters not in final_annotations) with NA for pipeline
  llm_raw_df$raw_pred[is.na(llm_raw_df$raw_pred)] <- NA_character_

  # ---- Build per-cell data frame ------------------------------------------
  test_cell_df <- data.frame(
    cell_id = colnames(seurat_test),
    cluster = as.character(seurat_test$Ground_Truth_Cluster),
    true_label = as.character(seurat_test$Ground_Truth_Celltype),
    stringsAsFactors = FALSE
  )

  gt_labels <- unique(as.character(seurat_test$Ground_Truth_Celltype))

  # ---- Run normalisation pipeline (post-processing) -----------------------
  cl_onto <- tryCatch(
    load_cell_ontology(),
    error = function(e) {
      warning("Could not load Cell Ontology: ", e$message)
      NULL
    }
  )

  norm_df <- tryCatch(
    run_normalisation_pipeline(llm_raw_df, gt_labels, cl_onto),
    error = function(e) {
      warning("Normalisation pipeline failed: ", e$message)
      NULL
    }
  )

  if (is.null(norm_df)) {
    return(default_return())
  }

  # ---- Diagnostic: per-cluster prediction mapping --------------------------
  # Join cluster ID back to GT celltype name so output is human-readable
  cat(sprintf("\n=== [mLLMCelltype] Prediction mapping (%d clusters) ===\n", nrow(norm_df)))
  if (exists("mapping_df")) {
    cluster_to_ct <- setNames(
      as.character(mapping_df$Ground_Truth_Celltype),
      as.character(mapping_df$Ground_Truth_Cluster)
    )
    norm_df_diag <- norm_df
    norm_df_diag$gt_celltype <- cluster_to_ct[as.character(norm_df_diag$cluster)]
    print(norm_df_diag[, c("cluster", "gt_celltype", "raw_pred", "normalised_pred", "strict", "lenient", "mapping_method")],
      row.names = FALSE
    )
  } else {
    print(norm_df[, c("cluster", "raw_pred", "normalised_pred", "strict", "lenient", "mapping_method")],
      row.names = FALSE
    )
  }
  cat(sprintf(
    "strict Unknown: %d/%d  |  lenient Unknown: %d/%d\n\n",
    sum(norm_df$strict == "Unknown"), nrow(norm_df),
    sum(norm_df$lenient == "Unknown"), nrow(norm_df)
  ))

  metrics_dual <- tryCatch(
    build_confusion_and_compute_metrics(test_cell_df, norm_df, gt_labels),
    error = function(e) {
      warning("build_confusion_and_compute_metrics failed: ", e$message)
      NULL
    }
  )

  # ---- Expand to cell level -----------------------------------------------
  cell_strict <- norm_df$strict[match(test_cell_df$cluster, norm_df$cluster)]
  cell_lenient <- norm_df$lenient[match(test_cell_df$cluster, norm_df$cluster)]
  cell_strict[is.na(cell_strict)] <- "Unknown"
  cell_lenient[is.na(cell_lenient)] <- "Unknown"

  wp_scores <- norm_df$wp_sim[match(test_cell_df$cluster, norm_df$cluster)]

  # Entropy-based confidence from consensus details
  consensus_details <- if (!is.null(consensus_results$initial_results)) {
    consensus_results$initial_results$consensus_results
  } else {
    warning("mLLMCelltype: initial_results missing from consensus output — entropy confidence unavailable")
    NULL
  }
  test_cluster_ids <- as.character(seurat_test$Ground_Truth_Cluster)
  cluster_confidence <- sapply(test_cluster_ids, function(x) {
    if (!is.null(consensus_details) && x %in% names(consensus_details) &&
      !is.null(consensus_details[[x]]$entropy) && length(models) > 1) {
      1 - consensus_details[[x]]$entropy / log(length(models))
    } else {
      0.5
    }
  })

  # Use wp_sim where available, fall back to entropy-based confidence
  confidence_scores <- ifelse(is.na(wp_scores), cluster_confidence, as.numeric(wp_scores))

  print(data.frame(
    cell_id    = test_cell_df$cell_id,
    strict     = cell_strict,
    lenient    = cell_lenient,
    true_label = test_cell_df$true_label
  ))

  # ---- Build cluster_mapping_df for CSV export ----------------------------
  cluster_gt_map <- unique(data.frame(
    cluster               = test_cell_df$cluster,
    ground_truth_celltype = test_cell_df$true_label,
    stringsAsFactors      = FALSE
  ))
  cluster_mapping_df <- merge(norm_df, cluster_gt_map, by = "cluster", all.x = TRUE)
  cluster_mapping_df <- cluster_mapping_df[, c(
    "cluster", "ground_truth_celltype",
    "raw_pred", "normalised_pred",
    "strict", "lenient",
    "mapping_method", "wp_sim"
  )]

  return(list(
    predictions = cell_strict,
    lenient_predictions = cell_lenient,
    true_labels = as.character(seurat_test$Ground_Truth_Celltype),
    confidence_scores = as.numeric(confidence_scores),
    cell_ids = colnames(seurat_test),
    strict_kappa = if (!is.null(metrics_dual)) metrics_dual$strict$kappa else NA_real_,
    lenient_kappa = if (!is.null(metrics_dual)) metrics_dual$lenient$kappa else NA_real_,
    kappa_gap = if (!is.null(metrics_dual)) metrics_dual$kappa_gap else NA_real_,
    hallucination_rate = if (!is.null(metrics_dual)) metrics_dual$hallucination_rate else NA_real_,
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb,
    cluster_mapping_df = cluster_mapping_df
  ))
}

# For backward compatibility
run_mLLMCelltype <- run_mLLMCelltype_function
