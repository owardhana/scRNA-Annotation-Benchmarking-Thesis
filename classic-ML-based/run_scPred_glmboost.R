# run_scPred.R
#################################################
# scPred Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scPred Cell Type Annotation Function
#' 
#' Purpose: Run scPred algorithm using machine learning-based classification
#' Inputs:
#'   - seurat_train: Training Seurat object (used to train ML model)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scPred's getFeatureSpace → trainModel → scPredict workflow
run_scPred_glmboost_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("scPred", quietly = TRUE)) {
    stop("scPred package not available. Please install it first.")
  }
  library(scPred)
  library(Seurat)
  library(caret)

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

  # Track peak memory usage
  runtime_secs <- NA
  peak_system_memory_mb <- NA

  # Run scPred workflow with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")

      start_time <- Sys.time()
    # 1. Create feature space from training data
    seurat_train <- getFeatureSpace(seurat_train, "Ground_Truth_Celltype")

    # 2. Train machine learning model
    seurat_train <- trainModel(seurat_train, model = "glmboost", allowParallel = TRUE)

    # 3. Predict on test data
    seurat_test <- scPredict(seurat_test, seurat_train)
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        # 1. Create feature space from training data
        seurat_train <- getFeatureSpace(seurat_train, "Ground_Truth_Celltype")

        # 2. Train machine learning model
        seurat_train <- trainModel(seurat_train, model = "glmboost", allowParallel = TRUE)

        # 3. Predict on test data
        seurat_test <- scPredict(seurat_test, seurat_train)
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }
  
  # Extract predictions and confidence scores
  cell_predictions <- seurat_test$scpred_prediction
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Extract confidence scores from probabilities stored in metadata
  prob_cols <- grep("scpred_", colnames(seurat_test@meta.data), value = TRUE)
  prob_cols <- prob_cols[prob_cols != "scpred_prediction"]
  
  if(length(prob_cols) > 0) {
    prob_matrix <- seurat_test@meta.data[, prob_cols, drop = FALSE]
    confidence_scores <- apply(prob_matrix, 1, max, na.rm = TRUE)
  } else {
    confidence_scores <- rep(1.0, length(cell_predictions))
    print("No confidence found")
    print(seurat_test@meta.data)
  }
  
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
    warning(paste("scPred (glmboost) error:", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scPred_glmboost <- run_scPred_glmboost_function
