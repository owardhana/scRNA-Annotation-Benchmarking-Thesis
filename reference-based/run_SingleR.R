# run_SingleR.R
#################################################
# SingleR Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' SingleR Cell Type Annotation Function
#' 
#' Purpose: Run SingleR algorithm using reference-based annotation with custom markers
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts to SingleCellExperiment, uses filtered markers, runs SingleR
run_SingleR_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("SingleR package not available. Please install it first.")
  }
  library(SingleR)
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package not available. Please install it first.")
  }
  library(SingleCellExperiment)
  
  # Convert markers from FindAllMarkers format to SingleR marker format
  convert_markers_to_singler <- function(marker_df) {
    # Filter for high quality markers first
    filtered_markers <- marker_df[marker_df$avg_log2FC >= 1.0 & 
                                    marker_df$p_val_adj < 0.05 & 
                                    marker_df$pct.1 >= 0.25, ]
    
    if(nrow(filtered_markers) == 0) {
      warning("No high-quality markers found")
      return(NULL)
    }
    
    # Create marker gene list by cluster
    marker_list <- split(filtered_markers$gene, filtered_markers$cluster)
    marker_list <- lapply(marker_list, unique)
    
    return(marker_list)
  }
  
  # Convert markers to SingleR format
  marker_list <- convert_markers_to_singler(markers)
  
  if(is.null(marker_list)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }
  
  # Ensure seurat_test has normalized data - use cross-compatible approach
  tryCatch({
    # Try to get normalized data - this will fail if not available
    test_data <- GetAssayData(seurat_test, assay = "RNA", layer = "data")
    if(length(test_data) == 0) {
      seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
    }
  }, error = function(e) {
    # Fallback for Seurat v4 or if layer approach fails
    tryCatch({
      test_data <- GetAssayData(seurat_test, assay = "RNA", layer = "data")
      if(length(test_data) == 0) {
        seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
      }
    }, error = function(e2) {
      # If both fail, normalize the data
      seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
    })
  })
  
  # Ensure seurat_train has normalized data - use cross-compatible approach
  tryCatch({
    # Try to get normalized data - this will fail if not available
    train_data <- GetAssayData(seurat_train, assay = "RNA", layer = "data")
    if(length(train_data) == 0) {
      seurat_train <- NormalizeData(seurat_train, verbose = FALSE)
    }
  }, error = function(e) {
    # Fallback for Seurat v4 or if layer approach fails
    tryCatch({
      train_data <- GetAssayData(seurat_train, assay = "RNA", layer = "data")
      if(length(train_data) == 0) {
        seurat_train <- NormalizeData(seurat_train, verbose = FALSE)
      }
    }, error = function(e2) {
      # If both fail, normalize the data
      seurat_train <- NormalizeData(seurat_train, verbose = FALSE)
    })
  })
  
  # Convert Seurat objects to SingleCellExperiment
  sce_train <- tryCatch({
    as.SingleCellExperiment(seurat_train)
  }, error = function(e) {
    warning(paste("Failed to convert training data to SingleCellExperiment:", e$message))
    return(NULL)
  })
  
  sce_test <- tryCatch({
    as.SingleCellExperiment(seurat_test)
  }, error = function(e) {
    warning(paste("Failed to convert test data to SingleCellExperiment:", e$message))
    return(NULL)
  })
  
  if(is.null(sce_train) || is.null(sce_test)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }

  # Validate data dimensions
  if(is.null(sce_train) || is.null(sce_test)) {
    warning("SingleCellExperiment conversion failed")
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }
  
  # Debug information
  cat(sprintf("Training data: %d genes x %d cells\n", nrow(sce_train), ncol(sce_train)))
  cat(sprintf("Test data: %d genes x %d cells\n", nrow(sce_test), ncol(sce_test)))
  cat(sprintf("Marker lists: %d clusters\n", length(marker_list)))

  # Track peak memory usage
  peak_memory_mb <- NA

  # Run SingleR with validated data and memory tracking
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("bench package not available for memory tracking")
    result <- tryCatch({
      SingleR(test = sce_test,
              ref = sce_train,
              labels = sce_train$Ground_Truth_Celltype,
              genes = marker_list,
              de.method = "wilcox")
    }, error = function(e) {
      warning(paste("SingleR execution failed:", e$message))
      return(NULL)
    })
  } else {
    library(bench)
    bench_result <- bench::mark(
      {
        result <- SingleR(test = sce_test,
                         ref = sce_train,
                         labels = sce_train$Ground_Truth_Celltype,
                         genes = marker_list,
                         de.method = "wilcox")
      },
      memory = TRUE,
      iterations = 1,
      check = FALSE
    )
    peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
  }
  
  if(is.null(result)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = peak_memory_mb
    ))
  }
  
  # Extract results - use result$labels (not result$label)
  if(!is.null(result$labels)) {
    cell_predictions <- result$labels
  } else {
    cell_predictions <- rep("Unknown", ncol(seurat_test))
  }
  
  # Use result$delta.max as confidence scores
  if(!is.null(result$delta.max)) {
    confidence_scores <- result$delta.max
  } else {
    confidence_scores <- rep(0, ncol(seurat_test))
  }
  
  # Handle NA predictions
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0
  
  # Validate results length
  if(length(cell_predictions) != ncol(seurat_test)) {
    warning(paste("Prediction length mismatch: expected", ncol(seurat_test), "got", length(cell_predictions)))
    cell_predictions <- rep("Unknown", ncol(seurat_test))
    confidence_scores <- rep(0, ncol(seurat_test))
  }
  
  # Get true labels from test set
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    peak_memory_mb = peak_memory_mb
  ))
}


# For backward compatibility
run_SingleR <- run_SingleR_function