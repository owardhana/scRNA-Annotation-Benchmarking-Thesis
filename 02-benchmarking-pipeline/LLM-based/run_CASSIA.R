# run_CASSIA.R
#################################################
# CASSIA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#
# Models:
#   annotation_model  : openai/gpt-5.4       (via OpenRouter)
#   score_model       : anthropic/claude-sonnet-4-6 (via OpenRouter)
#   annotationboost   : anthropic/claude-sonnet-4-6 (via OpenRouter)
#
# Post-processing: normalisation_pipeline.R maps cluster-level CASSIA
# annotations → strict/lenient GT labels (Wu-Palmer + 3-stage fallback).
#################################################

#' CASSIA Cell Type Annotation Function
#'
#' @param seurat_train Training Seurat object (used for marker context)
#' @param seurat_test  Test Seurat object to predict
#' @param markers      Marker genes dataframe from FindAllMarkers()
#'
#' @return List with predictions, true_labels, confidence_scores, cell_ids,
#'   plus lenient_predictions, strict_kappa, lenient_kappa, kappa_gap,
#'   hallucination_rate, runtime_secs, peak_system_memory_mb.
run_CASSIA_function <- function(seurat_train, seurat_test, markers) {

  # Use centralized prepare_markers() helper for uniform filtering
  source("benchmarking_helpers.R", local = TRUE)
  markers <- prepare_markers(markers)

  # Load normalisation pipeline
  source("LLM-based/normalisation_pipeline.R", local = TRUE)

  # Load required libraries with error checking
  required_packages <- c("CASSIA", "Seurat")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "not available. Please install it first."))
    }
  }

  library(CASSIA)
  library(Seurat)

  Idents(object = seurat_train) <- "Ground_Truth_Celltype"
  Idents(object = seurat_test)  <- "Ground_Truth_Celltype"

  # Set API keys
  openai_key     <- Sys.getenv("OPENAI_API_KEY")
  anthropic_key  <- Sys.getenv("ANTHROPIC_API_KEY")
  openrouter_key <- Sys.getenv("OPENROUTER_API_KEY")

  if (nchar(openai_key) > 0)     setLLMApiKey(openai_key,     provider = "openai")
  if (nchar(anthropic_key) > 0)  setLLMApiKey(anthropic_key,  provider = "anthropic")
  if (nchar(openrouter_key) > 0) {
    setLLMApiKey(openrouter_key, provider = "openrouter")
  } else {
    warning("OpenRouter API key not found - CASSIA requires valid API key")
  }

  # Create/use persistent CASSIA temp folder
  cassia_temp_dir <- file.path(getwd(), "misc", "CASSIA_temp")
  if (!dir.exists(cassia_temp_dir)) {
    dir.create(cassia_temp_dir, recursive = TRUE)
    warning("Created missing CASSIA_temp directory in misc folder")
  }
  if (!file.access(cassia_temp_dir, mode = 2) == 0) {
    stop("No write permission for CASSIA_temp directory: ", cassia_temp_dir)
  }

  old_wd <- getwd()
  setwd(cassia_temp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  output_name <- "CASSIA_results"

  # Default return
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

  cleanup <- function() {
    contents <- list.files(cassia_temp_dir, full.names = TRUE, all.files = FALSE)
    if (length(contents) > 0) unlink(contents, recursive = TRUE)
  }

  # Validate input
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in test data")
    cleanup()
    return(default_return())
  }

  runtime_secs         <- NA_real_
  peak_system_memory_mb <- NA_real_

  tryCatch({
    cat("Current working directory before CASSIA:", getwd(), "\n")
    cat("CASSIA_temp directory:", cassia_temp_dir, "\n")
    cat("Running CASSIA pipeline with output_name:", output_name, "\n")

    # ---- Run CASSIA pipeline (only this block is memory-tracked) ----------
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      runCASSIA_pipeline(
        output_file_name      = output_name,
        tissue                = "blood",
        species               = "human",
        marker                = markers,
        max_workers           = 4,
        annotation_model      = "openai/gpt-5.4",
        annotation_provider   = "openrouter",
        score_model           = "anthropic/claude-sonnet-4-6",
        score_provider        = "openrouter",
        annotationboost_model = "anthropic/claude-sonnet-4-6",
        annotationboost_provider = "openrouter",
        score_threshold       = 75,
        max_retries           = 1
      )
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        runCASSIA_pipeline(
          output_file_name      = output_name,
          tissue                = "blood",
          species               = "human",
          marker                = markers,
          max_workers           = 4,
          annotation_model      = "openai/gpt-5.4",
          annotation_provider   = "openrouter",
          score_model           = "anthropic/claude-sonnet-4-6",
          score_provider        = "openrouter",
          annotationboost_model = "anthropic/claude-sonnet-4-6",
          annotationboost_provider = "openrouter",
          score_threshold       = 75,
          max_retries           = 1
        )
      })
      runtime_secs         <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
      cat("Peak memory usage during CASSIA pipeline:", peak_system_memory_mb, "MB\n")
    }

    # ---- Locate CASSIA output CSV ------------------------------------------
    # CASSIA creates a timestamped subfolder inside CASSIA_temp; find the
    # newest one and read CASSIA_results_full.csv from it.
    subdirs <- list.dirs(cassia_temp_dir, recursive = FALSE, full.names = TRUE)
    if (length(subdirs) == 0) {
      warning("No CASSIA output subfolder found in CASSIA_temp.")
      cleanup()
      return(default_return())
    }
    newest_subdir <- subdirs[which.max(file.mtime(subdirs))]
    results_csv   <- file.path(newest_subdir, "CASSIA_results_full.csv")

    if (!file.exists(results_csv)) {
      warning("CASSIA_results_full.csv not found in ", newest_subdir)
      cleanup()
      return(default_return())
    }
    cat("Reading CASSIA results from:", results_csv, "\n")

    # ---- Read CSV and extract predictions per cluster --------------------
    cassia_df <- read.csv(results_csv, stringsAsFactors = FALSE)

    # True.Cell.Type corresponds to Ground_Truth_Cluster IDs.
    # Use them directly as the join key, matching the pattern in
    # GPTCelltype / mLLMCelltype_collab.
    llm_raw_df <- data.frame(
      cluster  = as.character(cassia_df$True.Cell.Type),
      raw_pred = as.character(cassia_df$Predicted.Main.Cell.Type),
      stringsAsFactors = FALSE
    )

    # ---- Build per-cell data frame ----------------------------------------
    # cluster = Ground_Truth_Cluster (join key), true_label = celltype name
    test_cell_df <- data.frame(
      cell_id    = colnames(seurat_test),
      cluster    = as.character(seurat_test$Ground_Truth_Cluster),
      true_label = as.character(seurat_test$Ground_Truth_Celltype),
      stringsAsFactors = FALSE
    )

    gt_labels <- unique(as.character(seurat_test$Ground_Truth_Celltype))

    # ---- Run normalisation pipeline (post-processing — not memory-tracked) -
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
      cleanup()
      return(default_return())
    }

    # ---- Diagnostic: per-cluster prediction mapping ------------------------
    cat(sprintf("\n=== [CASSIA] Prediction mapping (%d clusters) ===\n", nrow(norm_df)))
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

    # ---- Expand to cell level ---------------------------------------------
    cell_strict  <- norm_df$strict[match(test_cell_df$cluster, norm_df$cluster)]
    cell_lenient <- norm_df$lenient[match(test_cell_df$cluster, norm_df$cluster)]
    cell_strict[is.na(cell_strict)]   <- "Unknown"
    cell_lenient[is.na(cell_lenient)] <- "Unknown"

    wp_scores <- norm_df$wp_sim[match(test_cell_df$cluster, norm_df$cluster)]
    # Default to 1.0 where wp_sim is NA (exact/synonym/fuzzy matches)
    wp_scores_final <- ifelse(is.na(wp_scores), 1.0, wp_scores)

    cleanup()

    print(data.frame(
      cell_id    = test_cell_df$cell_id,
      strict     = cell_strict,
      lenient    = cell_lenient,
      true_label = test_cell_df$true_label
    ))

    # ---- Build cluster_mapping_df for CSV export ---------------------------
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
      confidence_scores    = as.numeric(wp_scores_final),
      cell_ids             = colnames(seurat_test),
      strict_kappa         = if (!is.null(metrics_dual)) metrics_dual$strict$kappa  else NA_real_,
      lenient_kappa        = if (!is.null(metrics_dual)) metrics_dual$lenient$kappa else NA_real_,
      kappa_gap            = if (!is.null(metrics_dual)) metrics_dual$kappa_gap     else NA_real_,
      hallucination_rate   = if (!is.null(metrics_dual)) metrics_dual$hallucination_rate else NA_real_,
      runtime_secs         = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb,
      cluster_mapping_df   = cluster_mapping_df
    ))

  }, error = function(e) {
    warning("CASSIA error: ", e$message)
    cleanup()
    return(default_return())
  })
}

# For backward compatibility
run_CASSIA <- run_CASSIA_function
