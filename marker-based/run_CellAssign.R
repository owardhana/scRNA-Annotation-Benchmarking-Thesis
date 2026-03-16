# run_CellAssign.R
#################################################
# CellAssign Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' Safe wrapper for cellassign that fixes TensorFlow shape conversion issues
#' 
#' This wrapper patches the Python-side reshape function to handle 
#' TensorShape objects that cause "Cannot convert a partially known TensorShape" errors
cellassign_safe <- function(...) {
  
  # Try to patch at Python level using reticulate
  patch_applied <- tryCatch({
    
    # Python code to monkey-patch tf.reshape
    py_code <- "
import tensorflow as tf

# Store original reshape
_original_reshape = tf.reshape

def _patched_reshape(tensor, shape, name=None):
    '''Patched reshape that handles partially known shapes'''
    try:
        # Try to convert shape to tensor if needed
        if hasattr(shape, '__iter__') and not isinstance(shape, (tf.Tensor, tf.Variable)):
            # Convert lists/tuples to tensor
            shape = tf.constant(shape, dtype=tf.int32)
    except:
        pass  # If conversion fails, use original shape
    
    return _original_reshape(tensor, shape, name=name)

# Apply the patch
tf.reshape = _patched_reshape
"
    
    # Execute the Python patch
    reticulate::py_run_string(py_code)
    
    TRUE
  }, error = function(e) {
    warning("Could not apply Python-level patch: ", e$message)
    FALSE
  })
  
  # Run cellassign
  result <- tryCatch({
    cellassign(...)
  }, error = function(e) {
    # If it still fails, provide helpful error message
    stop("CellAssign failed even with patch. Error: ", e$message,
         "\n\nTry these solutions:",
         "\n1. Reinstall cellassign from source with the fix",
         "\n2. Use TensorFlow 2.4.0 with tensorflow-probability 0.12.0",
         "\n3. Check GitHub issue #92 for manual file editing instructions")
  }, finally = {
    # Restore original reshape
    if (patch_applied) {
      tryCatch({
        reticulate::py_run_string("tf.reshape = _original_reshape")
      }, error = function(e) {
        # Silent fail on cleanup
      })
    }
  })
  
  return(result)
}

#' CellAssign Cell Type Annotation Function
#' 
#' Purpose: Run CellAssign algorithm using probabilistic assignment with TensorFlow
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict (raw counts required)
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts to SCE, creates marker list, uses raw counts + size factors
run_CellAssign_function <- function(seurat_train, seurat_test, markers) {
  
  # Simple fallback for any error
  unknown_result <- list(
    predictions = rep("Unknown", ncol(seurat_test)),
    true_labels = seurat_test$Ground_Truth_Celltype,
    confidence_scores = rep(0, ncol(seurat_test)),
    cell_ids = colnames(seurat_test)
  )
  
  tryCatch({
    # Load required packages
    library(SingleCellExperiment, quietly = TRUE)
    library(scuttle, quietly = TRUE)
    library(scran, quietly = TRUE)
    library(cellassign, quietly = TRUE)
    library(tensorflow, quietly = TRUE)
    
    # Convert Seurat to SCE
    sce <- as.SingleCellExperiment(seurat_test)
    
    # Create simple marker list (top 5 markers per cell type)
    marker_list <- list()
    cell_types <- unique(markers$cluster)
    
    for(ct in cell_types) {
      ct_markers <- markers[markers$cluster == ct & 
                              markers$avg_log2FC >= 1.0 & 
                              markers$p_val_adj < 0.05, "gene"]
      if(length(ct_markers) > 0) {
        marker_list[[ct]] <- head(ct_markers, 5)
      }
    }
    
    # Remove cell types with no markers
    marker_list <- marker_list[sapply(marker_list, length) > 0]
    
    if(length(marker_list) == 0) {
      return(unknown_result)
    }
    
    # Get all marker genes and subset SCE
    all_markers <- unique(unlist(marker_list))
    common_genes <- intersect(all_markers, rownames(sce))
    
    if(length(common_genes) < 5) {
      return(unknown_result)
    }
    
    sce_final <- sce[common_genes, ]
    
    # Simple size factors (library size normalization)
    s <- colSums(counts(sce_final))
    s <- s / median(s[s > 0])
    s[s <= 0] <- 1  # Replace zeros/negatives
    
    # Filter cells with sufficient expression
    valid_cells <- s > 0.1 & colSums(counts(sce_final)) >= 5
    
    if(sum(valid_cells) < 10) {
      return(unknown_result)
    }
    
    sce_final <- sce_final[, valid_cells]
    s <- s[valid_cells]

    # Run CellAssign with memory tracking (TensorFlow memory separate from R)
    peak_memory_mb <- NA
    if (!requireNamespace("bench", quietly = TRUE)) {
      warning("bench package not available for memory tracking")
      fit <- cellassign_safe(exprs_obj = sce_final,
                             marker_gene_info = marker_list,
                             s = s,
                             learning_rate = 1e-2,
                             shrinkage = TRUE,
                             verbose = FALSE)
    } else {
      library(bench)
      bench_result <- bench::mark(
        {
          fit <- cellassign_safe(exprs_obj = sce_final,
                                 marker_gene_info = marker_list,
                                 s = s,
                                 learning_rate = 1e-2,
                                 shrinkage = TRUE,
                                 verbose = FALSE)
        },
        memory = TRUE,
        iterations = 1,
        check = FALSE
      )
      # Note: This tracks R-side memory only; TensorFlow uses separate memory pool
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
    }

    # Extract results
    predictions_subset <- celltypes(fit)
    probabilities_subset <- cellprobs(fit)
    confidence_subset <- apply(probabilities_subset, 1, max)
    
    # Map back to all cells
    all_predictions <- rep("Unknown", ncol(seurat_test))
    names(all_predictions) <- colnames(seurat_test)
    
    all_confidence <- rep(0, ncol(seurat_test))
    names(all_confidence) <- colnames(seurat_test)
    
    # Map CellAssign results to valid cells
    valid_cell_names <- colnames(sce_final)
    all_predictions[valid_cell_names] <- as.character(predictions_subset)
    all_confidence[valid_cell_names] <- confidence_subset
    
    # Return standard format
    return(list(
      predictions = as.character(all_predictions),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = all_confidence,
      cell_ids = colnames(seurat_test),
      peak_memory_mb = peak_memory_mb
    ))
    
  }, error = function(e) {
    warning(paste("CellAssign failed:", e$message))
    return(unknown_result)
  })
}

# For backward compatibility
run_CellAssign <- run_CellAssign_function