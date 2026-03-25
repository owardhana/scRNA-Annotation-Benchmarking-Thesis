# run_SCINA.R
#################################################
# SCINA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' SCINA Cell Type Annotation Function
#' 
#' Purpose: Run SCINA algorithm on test data using training data signatures
#' Inputs:
#'    - seurat_train: Training Seurat object with cell type labels
#'    - seurat_test: Test Seurat object to predict
#'    - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts markers to SCINA signatures, runs SCINA, extracts results
run_SCINA_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("SCINA", quietly = TRUE)) {
    stop("SCINA package not available. Please install it first.")
  }
  library(SCINA)
  
  # Convert markers from FindAllMarkers format to SCINA signatures format
  convert_markers_to_scina <- function(marker_df, markers_per_type = 50) {
    # Validate input
    if(is.null(marker_df) || nrow(marker_df) == 0) {
      warning("marker_df is NULL or empty")
      return(NULL)
    }
    
    # Check for required columns
    required_cols <- c("cluster", "gene", "avg_log2FC")
    if(!all(required_cols %in% colnames(marker_df))) {
      warning(sprintf("marker_df missing required columns: %s",
                      paste(setdiff(required_cols, colnames(marker_df)), collapse = ", ")))
      return(NULL)
    }
    
    # Remove rows with NA cluster values (can't assign to a cell type)
    marker_df <- marker_df[!is.na(marker_df$cluster), ]
    if(nrow(marker_df) == 0) {
      warning("No valid markers after removing NA clusters")
      return(NULL)
    }
    
    signatures <- list()
    cell_types <- unique(marker_df$cluster)
    
    cat("\n=== Converting markers to SCINA format ===\n")
    cat(sprintf("Total markers from FindAllMarkers: %d\n", nrow(marker_df)))
    cat(sprintf("Cell types: %d\n", length(cell_types)))
    
    for(ct in cell_types) {
      # Get genes for this cell type
      ct_markers <- marker_df[marker_df$cluster == ct, ]
      
      # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
      filtered <- ct_markers[ct_markers$avg_log2FC >= 0.5 &
                               ct_markers$p_val_adj < 0.05 &
                               ct_markers$pct.1 >= 0.15, ]
      
      # Remove rows with NA in critical columns (keep this safeguard)
      filtered <- filtered[!is.na(filtered$avg_log2FC) & !is.na(filtered$gene), ]
      
      # If no valid markers for this cell type, skip it
      if(nrow(filtered) == 0) {
        warning(sprintf("No valid markers for cell type: %s", ct))
        next
      }
      
      # Sort by avg_log2FC and take top 50 (standardized)
      filtered <- filtered[order(filtered$avg_log2FC, decreasing = TRUE, na.last = TRUE), ]
      top_markers <- head(filtered$gene, 50)
      
      signatures[[as.character(ct)]] <- top_markers
      cat(sprintf("  %s: %d markers\n", ct, length(top_markers)))
    }
    
    return(signatures)
  }
  
  cat("\n=== SCINA: Semi-supervised Cell Type Annotation ===\n")
  
  # Extract expression data from test set
  # SCINA requires: genes (rows) x cells (columns) expression matrix
  exp_data <- tryCatch({
    as.matrix(GetAssayData(seurat_test, layer = "data"))
  }, error = function(e) {
    as.matrix(GetAssayData(seurat_test, slot = "data"))  # Fallback for Seurat v4
  })
  cat(sprintf("Expression data: %d genes x %d cells\n", nrow(exp_data), ncol(exp_data)))
  
  # Convert markers to SCINA signatures format
  signatures <- convert_markers_to_scina(markers, markers_per_type = 50)
  
  # Filter signatures to only include genes present in the expression data
  cat("\n=== Filtering signatures by gene overlap ===\n")
  signatures_filtered <- lapply(names(signatures), function(ct) {
    genes <- signatures[[ct]]
    overlap <- genes[genes %in% rownames(exp_data)]
    cat(sprintf("  %s: %d/%d genes overlap\n", ct, length(overlap), length(genes)))
    return(overlap)
  })
  names(signatures_filtered) <- names(signatures)
  
  # Remove signatures with too few genes (need at least 3 for SCINA to work)
  min_genes_required <- 3
  signatures_filtered <- signatures_filtered[sapply(signatures_filtered, length) >= min_genes_required]
  
  cat(sprintf("\nFinal signatures after filtering: %d cell types\n", length(signatures_filtered)))
  cat("Genes per signature:\n")
  for(ct in names(signatures_filtered)) {
    cat(sprintf("  %s: %d genes\n", ct, length(signatures_filtered[[ct]])))
  }
  
  if(is.null(signatures_filtered) || length(signatures_filtered) == 0) {
    warning("No valid signatures found for SCINA (all had < 3 genes)")
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test)
    ))
  }
  
  if(length(signatures_filtered) < 2) {
    warning("SCINA requires at least 2 cell types with valid signatures")
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test)
    ))
  }
  
  # Run SCINA with memory tracking
  cat("\n=== Running SCINA algorithm ===\n")
  cat("Parameters:\n")
  cat("  max_iter: 100\n")
  cat("  convergence_n: 10\n")
  cat("  convergence_rate: 0.99\n")
  cat("  rm_overlap: TRUE\n")
  cat("  allow_unknown: FALSE\n\n")
  
  runtime_secs <- NA
  peak_system_memory_mb <- NA
  
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
    start_time <- Sys.time()
    results <- tryCatch({
      SCINA(exp_data,
            signatures_filtered,
            max_iter = 100,
            convergence_n = 10,
            convergence_rate = 0.99,
            sensitivity_cutoff = 1,
            rm_overlap = TRUE,
            allow_unknown = FALSE
      )
    }, error = function(e) {
      warning(paste("SCINA failed:", e$message))
      cat("Error details:", e$message, "\n")
      return(NULL)
    })
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM({
      # Use <<- so 'results' is available outside the peakRAM environment
      results <<- SCINA(
        exp_data,
        signatures_filtered,
        max_iter = 100,
        convergence_n = 10,
        convergence_rate = 0.99,
        sensitivity_cutoff = 1,
        rm_overlap = TRUE,
        allow_unknown = FALSE
      )
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }
  
  if(is.null(results)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }
  
  # Extract results
  cat("\n=== SCINA completed successfully ===\n")
  cell_predictions <- results$cell_labels
  cell_probabilities <- results$probabilities
  
  # Get confidence scores (maximum probability for each cell)
  confidence_scores <- apply(cell_probabilities, 2, max)
  
  # Get true labels from test set
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Summary
  cat(sprintf("\nPrediction summary:\n"))
  cat(sprintf("  Total cells: %d\n", length(cell_predictions)))
  cat(sprintf("  Assigned: %d\n", sum(cell_predictions != "unknown")))
  cat(sprintf("  Unknown: %d\n", sum(cell_predictions == "unknown")))
  cat(sprintf("  Unique predicted types: %d\n",
              length(unique(cell_predictions[cell_predictions != "unknown"]))))
  cat(sprintf("  Mean confidence: %.3f\n", mean(confidence_scores)))
  cat(sprintf("  Median confidence: %.3f\n", median(confidence_scores)))
  
  cat("\nPredictions distribution:\n")
  print(table(cell_predictions))
  
  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
}

# For backward compatibility, keep the function available
run_SCINA <- run_SCINA_function