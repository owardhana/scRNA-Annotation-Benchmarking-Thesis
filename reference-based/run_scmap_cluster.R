# run_scmap_cluster.R
#################################################
# scmap Cluster Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scmap Cluster Cell Type Annotation Function
#' 
#' Purpose: Run scmap algorithm using reference-based cluster mapping
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scmapCluster for direct cluster-level cell type assignment
run_scmap_cluster_function <- function(seurat_train, seurat_test, markers) {
  
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
  Seurat_test_sce  <- as.SingleCellExperiment(seurat_test)
  
  # Add feature_symbol column required by scmap
  rowData(Seurat_train_sce)$feature_symbol <- make.unique(rownames(Seurat_train_sce))
  rowData(Seurat_test_sce)$feature_symbol  <- make.unique(rownames(Seurat_test_sce))
  
  # Ensure cluster labels come from metadata column "Ground_Truth_Celltype"
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
    stop("Training Seurat object must contain a 'Ground_Truth_Celltype' column in metadata")
  }
  colData(Seurat_train_sce)$cluster <- seurat_train$Ground_Truth_Celltype

  # Track peak memory usage
  peak_memory_mb <- NA

  # Run scmap cluster-level mapping with memory tracking
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("bench package not available for memory tracking")
    # Select features and create cluster index (explicitly using 'cluster' colData)
    Seurat_train_sce <- selectFeatures(Seurat_train_sce, suppress_plot = TRUE)
    Seurat_train_sce <- indexCluster(Seurat_train_sce, cluster_col = "cluster")

    # Run scmapCluster
    scmapCluster_results <- scmapCluster(
      projection = Seurat_test_sce,
      index_list = list(Seurat_train = metadata(Seurat_train_sce)$scmap_cluster_index)
    )
  } else {
    library(bench)
    bench_result <- bench::mark(
      {
        # Select features and create cluster index (explicitly using 'cluster' colData)
        Seurat_train_sce <- selectFeatures(Seurat_train_sce, suppress_plot = TRUE)
        Seurat_train_sce <- indexCluster(Seurat_train_sce, cluster_col = "cluster")

        # Run scmapCluster
        scmapCluster_results <- scmapCluster(
          projection = Seurat_test_sce,
          index_list = list(Seurat_train = metadata(Seurat_train_sce)$scmap_cluster_index)
        )
      },
      memory = TRUE,
      iterations = 1,
      check = FALSE
    )
    peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
  }
  
  # Extract direct cluster predictions (first column of matrix)
  predicted_labels <- scmapCluster_results$scmap_cluster_labs[, 1]

  # Add predictions back to test object
  colData(Seurat_test_sce)$predicted_celltype <- predicted_labels

  # Extract results
  cell_predictions <- predicted_labels
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Get confidence scores from scmapCluster similarities
  if ("scmap_cluster_siml" %in% names(scmapCluster_results)) {
    confidence_scores <- apply(scmapCluster_results$scmap_cluster_siml, 1, max, na.rm = TRUE)
  } else {
    confidence_scores <- rep(1.0, length(cell_predictions))
  }

  # Handle NA and "unassigned" predictions (scmap uses "unassigned" for low-confidence cells)
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  cell_predictions[cell_predictions == "unassigned"] <- "Unknown"
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
    warning(sprintf("scmap cluster failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scmap_cluster <- run_scmap_cluster_function