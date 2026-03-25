# run_ScInfeR.R
#################################################
# ScInfeR Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' ScInfeR Cell Type Annotation Function
#'
#' Purpose: Run ScInfeR algorithm using marker-based prediction with k-nearest neighbors
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict (will create clusters if needed)
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts markers to ScInfeR format (celltype, marker, weight),
#'            runs ScInfeR prediction with k-NN weighting
#' Note: ScInfeR uses 3-column marker format with weights (all set to 1)
run_ScInfeR_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("ScInfeR", quietly = TRUE)) {
    stop("ScInfeR package not available. Please install it with: devtools::install_github('xuyungang/ScInfeR')")
  }
  library(ScInfeR)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required. Please install it with: install.packages('dplyr')")
  }
  library(dplyr)

  # WORKAROUND: ScInfeR has a bug where it calls log_warning() but doesn't define it
  # Define it in global environment if it doesn't exist
  log_warning_added <- FALSE
  if (!exists("log_warning", envir = .GlobalEnv)) {
    log_warning <- function(...) {
      warning(paste(..., sep = " "))
    }
    assign("log_warning", log_warning, envir = .GlobalEnv)
    log_warning_added <- TRUE
  }

  # Ensure cleanup happens even if function exits early
  on.exit({
    if (log_warning_added && exists("log_warning", envir = .GlobalEnv)) {
      rm(log_warning, envir = .GlobalEnv)
    }
  })

  # Default return for errors
  default_return <- function() {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }

  # Validate inputs
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in metadata")
    return(default_return())
  }

  if (!is.data.frame(markers) || nrow(markers) == 0) {
    warning("Invalid or empty markers dataframe")
    return(default_return())
  }

  cat("\n=== ScInfeR: Marker-Based Annotation with k-NN ===\n")

  # Step 1: Ensure test data has PCA, UMAP, and clusters (ScInfeR requires UMAP)
  tryCatch({
    if (!"pca" %in% names(seurat_test@reductions)) {
      cat("Computing PCA for test data...\n")
      seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
      seurat_test <- FindVariableFeatures(seurat_test, verbose = FALSE)
      seurat_test <- ScaleData(seurat_test, verbose = FALSE)
      seurat_test <- RunPCA(seurat_test, verbose = FALSE)
    }
    if (!"umap" %in% names(seurat_test@reductions)) {
      cat("Computing UMAP for test data...\n")
      seurat_test <- RunUMAP(seurat_test, dims = 1:min(30, ncol(seurat_test[["pca"]])), verbose = FALSE)
    }
    if (!"seurat_clusters" %in% colnames(seurat_test@meta.data)) {
      cat("Creating clusters for test data...\n")
      seurat_test <- FindNeighbors(seurat_test, verbose = FALSE)
      seurat_test <- FindClusters(seurat_test, resolution = 0.5, verbose = FALSE)
    }
  }, error = function(e) {
    warning(paste("Failed to prepare test data:", e$message))
    return(default_return())
  })

  n_test_clusters <- length(unique(seurat_test$seurat_clusters))
  cat(sprintf("Test data: %d clusters across %d cells\n",
              n_test_clusters, ncol(seurat_test)))

  # Step 2: Convert markers to ScInfeR format (3 columns: celltype, marker, weight)
  cat("\nConverting markers to ScInfeR format (celltype, marker, weight)...\n")

  # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
  filtered_markers <- markers[markers$avg_log2FC >= 0.5 &
                              markers$p_val_adj < 0.05 &
                              markers$pct.1 >= 0.15, ]

  if (nrow(filtered_markers) == 0) {
    warning("No high-quality markers found")
    return(default_return())
  }

  cat(sprintf("Filtered markers: %d genes across %d cell types\n",
              nrow(filtered_markers),
              length(unique(filtered_markers$cluster))))

  # Take top 50 markers per cell type (standardized)
  markers_per_type <- 50
  top_markers <- filtered_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = markers_per_type) %>%
    dplyr::ungroup()

  cat(sprintf("Using top %d markers per cell type (standardized)\n", markers_per_type))
  cat("Markers per cell type:\n")
  print(table(top_markers$cluster))

  # Convert to ScInfeR 3-column format
  scinfer_markers <- data.frame(
    celltype = as.character(top_markers$cluster),
    marker = top_markers$gene,
    weight = 1,  # All weights set to 1
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nScInfeR marker dataframe created:\n"))
  cat(sprintf("  Rows: %d (marker entries)\n", nrow(scinfer_markers)))
  cat(sprintf("  Columns: 3 (celltype, marker, weight)\n"))
  cat(sprintf("  Cell types: %d\n", length(unique(scinfer_markers$celltype))))
  cat(sprintf("  All weights set to: 1\n"))

  cat("\nSample markers:\n")
  print(head(scinfer_markers, 10))

  # Step 3: Run ScInfeR prediction
  cat("\nRunning ScInfeR prediction...\n")
  cat("Parameters:\n")
  cat("  own_weightage: 0.5\n")
  cat("  n_neighbor: 10\n")
  cat("  assay: RNA\n")
  cat("  slot: counts (raw counts)\n")
  cat("  Note: Using 'counts' slot as ScInfeR may have issues with 'data' slot\n\n")

  # Increase C stack limit to handle ScInfeR's recursion
  old_cstack <- Cstack_info()["size"]
  tryCatch({
    # Try to increase stack size (may not work on all systems)
    options(expressions = 500000)
  }, error = function(e) {
    cat("Note: Could not increase expression limit\n")
  })

  runtime_secs <- NA
  peak_system_memory_mb <- NA
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
    scinfer_result <- tryCatch({
      predict_celltype_scRNA_seurat(
        s_object = seurat_test,
        group_annt = seurat_test$seurat_clusters,
        ct_marker_df = scinfer_markers,
        subtype_present = FALSE,
        subtype_info = FALSE,
        assay_name = "RNA",
        slot_name = "counts",
        own_weightage = 0.5,
        n_neighbor = 10
      )
    }, error = function(e) {
      warning(paste("ScInfeR prediction failed:", e$message))
      cat("Error details:", e$message, "\n")
      return(NULL)
    })
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        scinfer_result <- predict_celltype_scRNA_seurat(
          s_object = seurat_test,
          group_annt = seurat_test$seurat_clusters,
          ct_marker_df = scinfer_markers,
          subtype_present = FALSE,
          subtype_info = FALSE,
          assay_name = "RNA",
          slot_name = "counts",
          own_weightage = 0.5,
          n_neighbor = 10
        )
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  if (is.null(scinfer_result)) {
    return(default_return())
  }

  # Step 4: Extract predictions
  if (!"celltype" %in% names(scinfer_result)) {
    warning("ScInfeR did not return 'celltype' column in results")
    return(default_return())
  }

  predictions <- as.character(scinfer_result$celltype)
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Validate predictions length
  if (length(predictions) != ncol(seurat_test)) {
    warning(sprintf("Prediction count mismatch: got %d, expected %d",
                    length(predictions), ncol(seurat_test)))
    return(default_return())
  }

  # Step 5: Get confidence scores
  # ScInfeR doesn't return confidence scores explicitly
  # Set all to 1.0 for successful predictions, 0.0 for Unknown
  confidence_scores <- ifelse(is.na(predictions) | predictions == "Unknown", 0, 1)

  # Handle NA predictions
  predictions[is.na(predictions)] <- "Unknown"

  # Summary
  cat(sprintf("\nPrediction summary:\n"))
  cat(sprintf("  Total cells: %d\n", length(predictions)))
  cat(sprintf("  Assigned: %d\n", sum(predictions != "Unknown")))
  cat(sprintf("  Unknown: %d\n", sum(predictions == "Unknown")))
  cat(sprintf("  Unique predicted types: %d\n",
              length(unique(predictions[predictions != "Unknown"]))))

  cat("\nPredictions distribution:\n")
  print(table(predictions))

  # Return standardized format
  return(list(
    predictions = as.character(predictions),
    true_labels = true_labels,
    confidence_scores = as.numeric(confidence_scores),
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
}

# For backward compatibility
run_ScInfeR <- run_ScInfeR_function
