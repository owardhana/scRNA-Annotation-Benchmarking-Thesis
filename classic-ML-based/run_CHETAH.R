# run_CHETAH.R
#################################################
# CHETAH Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' CHETAH Cell Type Annotation Function
#' 
#' Purpose: Run CHETAH algorithm using reference-based classification
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses CHETAH's CHETAHclassifier with SingleCellExperiment objects
run_CHETAH_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required libraries
  if (!requireNamespace("CHETAH", quietly = TRUE)) {
    stop("CHETAH package not available. Please install it first.")
  }
  library(CHETAH)
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package not available. Please install it first.")
  }
  library(SingleCellExperiment)

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

  # 1. Convert Seurat objects to SingleCellExperiment
  sce_train <- as.SingleCellExperiment(seurat_train)
  sce_test <- as.SingleCellExperiment(seurat_test)

  # Add cell type data to colData of reference
  # CHETAH requires that "Unassigned" not be used in cell type names
  # Replace "Unassigned" with "unassigned" (lowercase) to comply
  cell_types <- seurat_train$Ground_Truth_Celltype
  cell_types[cell_types == "Unassigned"] <- "unassigned"
  colData(sce_train)$cell_type <- cell_types

  # Track peak memory usage
  runtime_secs <- NA
  peak_system_memory_mb <- NA

  # 2. Run CHETAH classifier with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
    sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce_train, ref_ct = "cell_type")
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce_train, ref_ct = "cell_type")
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  # 3. Extract predictions
  cell_predictions <- sce_test$celltype_CHETAH
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Convert predictions back to "Unassigned" (uppercase) to match ground truth
  cell_predictions[cell_predictions == "unassigned"] <- "Unassigned"

  # 4. Extract confidence scores from nested DataFrame
  if(!is.null(sce_test$int_colData$CHETAH) && "conf_scores" %in% names(sce_test$int_colData$CHETAH)) {
    confidence_scores <- sce_test$int_colData$CHETAH$conf_scores
    # If conf_scores is a list, extract appropriate values
    if(is.list(confidence_scores)) {
      confidence_scores <- sapply(confidence_scores, function(x) max(x, na.rm = TRUE))
    }
  } else {
    confidence_scores <- rep(1.0, length(cell_predictions))
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
    warning(paste("CHETAH error:", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_CHETAH <- run_CHETAH_function