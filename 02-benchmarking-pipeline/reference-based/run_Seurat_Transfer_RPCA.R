# run_Seurat_Transfer_RPCA.R
#################################################
# Seurat Transfer Learning RPCA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' Seurat Transfer Learning RPCA Cell Type Annotation Function
#' 
#' Purpose: Run Seurat transfer learning using RPCA reduction
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses FindTransferAnchors with RPCA reduction and TransferData
run_Seurat_Transfer_RPCA_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package not available. Please install it first.")
  }
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
    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Run Seurat transfer learning with RPCA and memory tracking
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      # Find transfer anchors using RPCA
      anchors <- FindTransferAnchors(
        reference = seurat_train,
        query = seurat_test,
        dims = 1:30,
        reduction = "rpca"
      )

      # Transfer data
      predictions <- TransferData(
        anchorset = anchors,
        refdata = seurat_train$Ground_Truth_Celltype,
        dims = 1:30,
        weight.reduction = "rpca.ref"
      )

      # Add predictions to seurat_test metadata
      seurat_test <- AddMetaData(seurat_test, metadata = predictions)
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        # Find transfer anchors using RPCA
        anchors <- FindTransferAnchors(
          reference = seurat_train,
          query = seurat_test,
          dims = 1:30,
          reduction = "rpca"
        )

        # Transfer data
        predictions <- TransferData(
          anchorset = anchors,
          refdata = seurat_train$Ground_Truth_Celltype,
          dims = 1:30,
          weight.reduction = "rpca.ref"
        )

        # Add predictions to seurat_test metadata
        seurat_test <- AddMetaData(seurat_test, metadata = predictions)
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    # Extract results - seurat_test$predicted.id vs seurat_test$Ground_Truth_Celltype
    cell_predictions <- seurat_test$predicted.id
    true_labels <- seurat_test$Ground_Truth_Celltype
    confidence_scores <- seurat_test@meta.data$prediction.score.max

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
    warning(sprintf("Seurat Transfer RPCA failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_Seurat_Transfer_RPCA <- run_Seurat_Transfer_RPCA_function