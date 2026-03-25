# run_scClassify.R
#################################################
# scClassify Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scClassify Cell Type Annotation Function
#' 
#' Purpose: Run scClassify algorithm using ensemble machine learning classification
#' Inputs:
#'   - seurat_train: Training Seurat object (used to train classifier)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scClassify with HOPACH tree, WKNN algorithm, limma features, cosine similarity
run_scClassify_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("scClassify", quietly = TRUE)) {
    stop("scClassify package not available. Please install it first.")
  }
  library(scClassify)
  library(Seurat)

  # Default return for error handling
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

  tryCatch({

  # Extract log-transformed matrices (genes x cells) - cross-compatible approach
  seurat_train_mtx <- GetAssayData(seurat_train, assay = "RNA", layer = "data")

  seurat_test_mtx <- GetAssayData(seurat_test, assay = "RNA", layer = "data")

  # Track peak memory usage
  runtime_secs <- NA
  peak_system_memory_mb <- NA

  # Run scClassify with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
    scClassify_res <- scClassify(
      exprsMat_train = seurat_train_mtx,
      cellTypes_train = seurat_train$Ground_Truth_Celltype,
      exprsMat_test = list(test = seurat_test_mtx),
      tree = "HOPACH",
      algorithm = "WKNN",
      selectFeatures = c("limma"),
      similarity = c("cosine"),
      returnList = FALSE,
      verbose = FALSE
    )
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        scClassify_res <- scClassify(
          exprsMat_train = seurat_train_mtx,
          cellTypes_train = seurat_train$Ground_Truth_Celltype,
          exprsMat_test = list(test = seurat_test_mtx),
          tree = "HOPACH",
          algorithm = "WKNN",
          selectFeatures = c("limma"),
          similarity = c("cosine"),
          returnList = FALSE,
          verbose = FALSE
        )
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  # Extract predictions from nested result structure
  cell_predictions <- scClassify_res$testRes$test$cosine_WKNN_limma$predRes
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Create confidence scores (scClassify doesn't provide explicit confidence)
  confidence_scores <- rep(1.0, length(cell_predictions))

  # Handle NA values
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))

  }, error = function(e) {
    warning(paste("scClassify error:", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scClassify <- run_scClassify_function