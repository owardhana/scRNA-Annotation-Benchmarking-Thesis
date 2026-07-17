# run_GPTCelltype.R
#
# GPTCelltype cell type annotation using the normalisation pipeline.
# Model: gpt-5.4
#
# Predictions are mapped cluster-level (LLM text → GT label) via
# normalisation_pipeline.R, yielding both STRICT and LENIENT labels.
# The primary `predictions` field contains strict labels for framework
# compatibility; lenient_predictions and kappa statistics are additional fields.

run_GPTCelltype_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries with error checking
  required_packages <- c("GPTCelltype", "Seurat")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "not available. Please install it first."))
    }
  }
  library(GPTCelltype)
  library(Seurat)

  # Use centralized prepare_markers() helper for uniform filtering
  source("benchmarking_helpers.R", local = TRUE)
  #markers <- prepare_markers(markers)

  # Load normalisation pipeline
  source("LLM-based/normalisation_pipeline.R", local = TRUE)

  # Default return function for error handling
  default_return <- function() {
    n <- ncol(seurat_test)
    return(list(
      predictions          = rep("Unknown", n),
      lenient_predictions  = rep("Unknown", n),
      true_labels          = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores    = rep(0, n),
      cell_ids             = colnames(seurat_test),
      strict_kappa         = NA_real_,
      lenient_kappa        = NA_real_,
      kappa_gap            = NA_real_,
      hallucination_rate   = NA_real_,
      runtime_secs         = NA_real_,
      peak_system_memory_mb = NA_real_
    ))
  }

  # Validate API key before attempting LLM call
  openai_key <- Sys.getenv("OPENAI_API_KEY")
  if (nchar(openai_key) == 0) {
    warning("OPENAI_API_KEY not set — GPTCelltype requires a valid key")
    return(default_return())
  }

  # Track peak memory usage
  runtime_secs         <- NA_real_
  peak_system_memory_mb <- NA_real_

  # Run GPTCelltype — wrap only the API call in peakRAM
  res <- NULL
  print("starting predictions")
  start_time <- Sys.time()
  markers <- FindAllMarkers(seurat_train)
  res <- gptcelltype(markers, model = "gpt-4")
  runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  print("finished predictions")
  if (is.null(res)) return(default_return())

  # ---- Build cluster-level raw prediction data frame ----------------------
  # gptcelltype() returns a named character vector:
  #   names  = marker cluster names (= Ground_Truth_Celltype labels from FindAllMarkers)
  #   values = GPT free-text cell type names
  unique_clusters <- names(res)
  llm_raw_df <- data.frame(
    cluster  = unique_clusters,
    raw_pred = as.character(res[unique_clusters]),
    stringsAsFactors = FALSE
  )

  # ---- Build per-cell data frame ------------------------------------------
  # norm_df$cluster is keyed by whatever markers$cluster contained:
  #   CV mode          → Ground_Truth_Celltype names
  #   Full-dataset mode → Ground_Truth_Cluster IDs
  # test_cell_df$cluster must use the same column so match() succeeds.
  marker_cluster_ids <- unique(as.character(markers$cluster))
  gt_celltype_vals   <- unique(as.character(seurat_test$Ground_Truth_Celltype))
  use_cluster_col    <- !all(marker_cluster_ids %in% gt_celltype_vals) &&
                        "Ground_Truth_Cluster" %in% colnames(seurat_test@meta.data)
  test_cluster_col   <- if (use_cluster_col) as.character(seurat_test$Ground_Truth_Cluster)
                        else                  as.character(seurat_test$Ground_Truth_Celltype)

  test_cell_df <- data.frame(
    cell_id    = colnames(seurat_test),
    cluster    = test_cluster_col,
    true_label = as.character(seurat_test$Ground_Truth_Celltype),
    stringsAsFactors = FALSE
  )

  gt_labels <- unique(as.character(seurat_test$Ground_Truth_Celltype))

  # ---- Run normalisation pipeline -----------------------------------------
  cl_onto <- tryCatch(
    load_cell_ontology(),
    error = function(e) {
      warning("Could not load Cell Ontology: ", e$message,
              " — Wu-Palmer stage will be skipped.")
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

  if (is.null(norm_df)) return(default_return())

  # ---- Diagnostic: per-cluster prediction mapping --------------------------
  cat(sprintf("\n=== [GPTCelltype] Prediction mapping (%d clusters) ===\n", nrow(norm_df)))
  print(norm_df[, c("cluster", "raw_pred", "normalised_pred", "strict", "lenient", "mapping_method")],
        row.names = FALSE)
  cat(sprintf("strict Unknown: %d/%d  |  lenient Unknown: %d/%d\n\n",
              sum(norm_df$strict  == "Unknown"), nrow(norm_df),
              sum(norm_df$lenient == "Unknown"), nrow(norm_df)))

  metrics_dual <- tryCatch(
    build_confusion_and_compute_metrics(test_cell_df, norm_df, gt_labels),
    error = function(e) {
      warning("build_confusion_and_compute_metrics failed: ", e$message)
      NULL
    }
  )

  # ---- Expand cluster-level labels to cell level --------------------------
  cell_strict  <- norm_df$strict[match(test_cell_df$cluster, norm_df$cluster)]
  cell_lenient <- norm_df$lenient[match(test_cell_df$cluster, norm_df$cluster)]
  cell_strict[is.na(cell_strict)]   <- "Unknown"
  cell_lenient[is.na(cell_lenient)] <- "Unknown"

  wp_scores <- norm_df$wp_sim[match(test_cell_df$cluster, norm_df$cluster)]
  wp_scores[is.na(wp_scores)] <- 0

  print(data.frame(
    cell_id     = test_cell_df$cell_id,
    strict      = cell_strict,
    lenient     = cell_lenient,
    true_label  = test_cell_df$true_label
  ))

  # ---- Build cluster_mapping_df for CSV export ----------------------------
  cluster_gt_map <- unique(data.frame(
    cluster               = test_cell_df$cluster,
    ground_truth_celltype = test_cell_df$true_label,
    stringsAsFactors      = FALSE
  ))
  cluster_mapping_df <- merge(norm_df, cluster_gt_map, by = "cluster", all.x = TRUE)
  cluster_mapping_df <- cluster_mapping_df[, c("cluster", "ground_truth_celltype",
                                                "raw_pred", "normalised_pred",
                                                "strict", "lenient",
                                                "mapping_method", "wp_sim")]

  return(list(
    predictions          = cell_strict,
    lenient_predictions  = cell_lenient,
    true_labels          = as.character(seurat_test$Ground_Truth_Celltype),
    confidence_scores    = as.numeric(wp_scores),
    cell_ids             = colnames(seurat_test),
    strict_kappa         = if (!is.null(metrics_dual)) metrics_dual$strict$kappa  else NA_real_,
    lenient_kappa        = if (!is.null(metrics_dual)) metrics_dual$lenient$kappa else NA_real_,
    kappa_gap            = if (!is.null(metrics_dual)) metrics_dual$kappa_gap     else NA_real_,
    hallucination_rate   = if (!is.null(metrics_dual)) metrics_dual$hallucination_rate else NA_real_,
    runtime_secs         = runtime_secs,
    peak_system_memory_mb = NA,
    cluster_mapping_df   = cluster_mapping_df
  ))
}

run_GPTCelltype <- run_GPTCelltype_function
