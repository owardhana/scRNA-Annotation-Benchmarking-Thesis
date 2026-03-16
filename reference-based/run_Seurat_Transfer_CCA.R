# run_Seurat_Transfer_CCA.R
#################################################
# Seurat Transfer Learning CCA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' Seurat Transfer Learning CCA Cell Type Annotation Function
#' 
#' Purpose: Run Seurat transfer learning using CCA reduction
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses FindTransferAnchors with CCA reduction and TransferData
run_Seurat_Transfer_CCA_function <- function(seurat_train, seurat_test, markers) {
  
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
      peak_memory_mb = NA
    ))
  }

  tryCatch({
    # Track peak memory usage
    peak_memory_mb <- NA

    # Run Seurat transfer learning with CCA and memory tracking
    if (!requireNamespace("bench", quietly = TRUE)) {
      warning("bench package not available for memory tracking")
      # Find transfer anchors using CCA
      anchors <- FindTransferAnchors(
        reference = seurat_train,
        query = seurat_test,
        dims = 1:30,
        reduction = "cca"
      )

      # Transfer data
      predictions <- TransferData(
        anchorset = anchors,
        refdata = seurat_train$Ground_Truth_Celltype,
        dims = 1:30,
        weight.reduction = "cca"
      )

      # Add predictions to seurat_test metadata
      seurat_test <- AddMetaData(seurat_test, metadata = predictions)
    } else {
      library(bench)
      bench_result <- bench::mark(
        {
          # Find transfer anchors using CCA
          anchors <- FindTransferAnchors(
            reference = seurat_train,
            query = seurat_test,
            dims = 1:30,
            reduction = "cca"
          )

          # Transfer data
          predictions <- TransferData(
            anchorset = anchors,
            refdata = seurat_train$Ground_Truth_Celltype,
            dims = 1:30,
            weight.reduction = "cca"
          )

          # Add predictions to seurat_test metadata
          seurat_test <- AddMetaData(seurat_test, metadata = predictions)
        },
        memory = TRUE,
        iterations = 1,
        check = FALSE
      )
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
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
      peak_memory_mb = peak_memory_mb
    ))
  }, error = function(e) {
    warning(sprintf("Seurat Transfer CCA failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_Seurat_Transfer_CCA <- run_Seurat_Transfer_CCA_function