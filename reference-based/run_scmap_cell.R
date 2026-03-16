# run_scmap.R
#################################################
# scmap Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scmap Cell Type Annotation Function
#' 
#' Purpose: Run scmap algorithm using reference-based cell mapping
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scmapCell with majority vote for cell type assignment
run_scmap_cell_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required libraries
  if (!requireNamespace("scmap", quietly = TRUE)) {
    stop("scmap package not available. Please install it first.")
  }
  library(scmap)
  
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
      peak_memory_mb = NA
    ))
  }

  tryCatch({
  # Convert to SingleCellExperiment
  Seurat_train_sce <- as.SingleCellExperiment(seurat_train)
  Seurat_test_sce <- as.SingleCellExperiment(seurat_test)
  
  # Add feature_symbol column required by scmap
  rowData(Seurat_train_sce)$feature_symbol <- rownames(Seurat_train_sce)
  rowData(Seurat_test_sce)$feature_symbol  <- rownames(Seurat_test_sce)
  
  # Make sure rownames are unique (scmap requires unique feature symbols)
  rowData(Seurat_train_sce)$feature_symbol <- make.unique(rowData(Seurat_train_sce)$feature_symbol)
  rowData(Seurat_test_sce)$feature_symbol  <- make.unique(rowData(Seurat_test_sce)$feature_symbol)

  # Track peak memory usage
  peak_memory_mb <- NA

  # Run scmap cell-level mapping with memory tracking
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("bench package not available for memory tracking")
    # Select features and create index
    Seurat_train_sce <- selectFeatures(Seurat_train_sce, suppress_plot = TRUE)
    Seurat_train_sce <- indexCell(Seurat_train_sce)

    # Run scmapCell
    scmapCell_results <- scmapCell(
      projection = Seurat_test_sce,
      index_list = list(Seurat_train = metadata(Seurat_train_sce)$scmap_cell_index)
    )
  } else {
    library(bench)
    bench_result <- bench::mark(
      {
        # Select features and create index
        Seurat_train_sce <- selectFeatures(Seurat_train_sce, suppress_plot = TRUE)
        Seurat_train_sce <- indexCell(Seurat_train_sce)

        # Run scmapCell
        scmapCell_results <- scmapCell(
          projection = Seurat_test_sce,
          index_list = list(Seurat_train = metadata(Seurat_train_sce)$scmap_cell_index)
        )
      },
      memory = TRUE,
      iterations = 1,
      check = FALSE
    )
    peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
  }

  # Get majority-vote labels
  predicted_labels <- apply(scmapCell_results$Seurat_train$cells, 2, function(x) {
    ref_labels <- colData(Seurat_train_sce)$Ground_Truth_Celltype[x]
    if (length(ref_labels) == 0) return(NA)
    names(sort(table(ref_labels), decreasing = TRUE))[1]
  })

  # Add predictions back to test object
  colData(Seurat_test_sce)$predicted_celltype <- predicted_labels

  # Extract results
  cell_predictions <- predicted_labels
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Get confidence scores from scmapCell similarities
  if("similarities" %in% names(scmapCell_results$Seurat_train)) {
    confidence_scores <- apply(scmapCell_results$Seurat_train$similarities, 2, max, na.rm = TRUE)
  } else {
    confidence_scores <- rep(1.0, length(cell_predictions))
  }
  
  # Handle NA predictions and confidence scores
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    peak_memory_mb = peak_memory_mb
  ))
  }, error = function(e) {
    warning(sprintf("scmap cell failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scmap <- run_scmap_cell_function