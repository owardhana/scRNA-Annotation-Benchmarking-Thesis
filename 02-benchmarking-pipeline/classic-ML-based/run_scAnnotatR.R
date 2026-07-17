# run_scAnnotatR.R

#' scAnnotatR Function (Seurat Input -> SCE Execution)
#' 
#' Purpose: Accepts Seurat objects, converts them to SCE, and runs scAnnotatR 
#' using its native SingleCellExperiment support.
#' 
#' Inputs:
#'    - seurat_train: Training Seurat object
#'    - seurat_test: Test Seurat object
#'    - markers: Marker genes dataframe
#' Outputs: Standardized list of predictions
run_scAnnotatR_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required libraries
  requireNamespace("scAnnotatR", quietly = TRUE)
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("SingleCellExperiment", quietly = TRUE)
  
  library(scAnnotatR)
  library(Seurat)
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
    # -------------------------------------------------------
    # 1. CONVERSION: Seurat -> SingleCellExperiment
    # -------------------------------------------------------
    # We convert the input Seurat objects to SCE format.
    # This preserves counts, logcounts, and metadata (colData).
    sce_train <- as.SingleCellExperiment(seurat_train)
    sce_test  <- as.SingleCellExperiment(seurat_test)
    
    # Ensure assay names are standard (scAnnotatR often looks for 'counts')
    # Seurat conversion usually creates 'counts' and 'logcounts'.
    # We verify 'counts' exists; if not, we try to alias it.
    if (!"counts" %in% assayNames(sce_train)) {
      if ("RNA" %in% assayNames(sce_train)) {
        assay(sce_train, "counts") <- assay(sce_train, "RNA")
      }
    }
    
    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA
    
    # Get unique cell types from the converted SCE object
    unique_celltypes <- unique(colData(sce_train)$Ground_Truth_Celltype)
    
    # Helper to run the train+predict pipeline
    run_scAnnotatR_pipeline <- function() {
      # -------------------------------------------------------
      # 2. TRAINING (Using SCE)
      # -------------------------------------------------------
      classifiers <- list()
      
      for (ct in unique_celltypes) {
        tryCatch({
          ct_df <- markers[markers$cluster == ct & markers$avg_log2FC > 0, ]
          ct_df <- ct_df[order(ct_df$avg_log2FC, decreasing = TRUE), ]
          ct_markers <- head(ct_df$gene, 20)
          available_genes <- intersect(ct_markers, rownames(sce_train))

          if(length(available_genes) >= 5) {
            classifier <- train_classifier(
              train_obj = sce_train,
              assay = "counts",
              cell_type = ct,
              marker_genes = available_genes,
              tag_slot = "Ground_Truth_Celltype"
            )
            classifiers[[ct]] <- classifier
          }
        }, error = function(e) {})
      }
      
      # -------------------------------------------------------
      # 3. PREDICTION (Using SCE)
      # -------------------------------------------------------
      if(length(classifiers) == 0) {
        return(NULL)
      }
      
      classify_cells(
        classify_obj = sce_test,
        classifiers = classifiers,
        assay = "counts"
      )
    }
    
    # Run scAnnotatR workflow with memory tracking
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      sce_test_classified <- run_scAnnotatR_pipeline()
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        sce_test_classified <- run_scAnnotatR_pipeline()
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }
    
    if(is.null(sce_test_classified)) {
      return(list(
        predictions = rep("Unknown", ncol(sce_test)),
        true_labels = as.character(colData(sce_test)$Ground_Truth_Celltype),
        confidence_scores = rep(0, ncol(sce_test)),
        cell_ids = colnames(sce_test),
        runtime_secs = runtime_secs,
        peak_system_memory_mb = peak_system_memory_mb
      ))
    }
    
    # -------------------------------------------------------
    # 4. EXTRACT RESULTS (From SCE colData)
    # -------------------------------------------------------
    results_meta <- colData(sce_test_classified)
    
    # Extract Predictions
    # scAnnotatR typically populates 'predicted_cell_type' or 'most_probable_cell_type'
    if("predicted_cell_type" %in% colnames(results_meta)) {
      preds <- results_meta$predicted_cell_type
    } else if("most_probable_cell_type" %in% colnames(results_meta)) {
      preds <- results_meta$most_probable_cell_type
    } else {
      preds <- rep("Unknown", ncol(sce_test))
    }
    
    # Extract Confidence Scores (Max probability across trained types)
    prob_cols <- grep("_p$", colnames(results_meta), value = TRUE)
    
    if(length(prob_cols) > 0) {
      prob_matrix <- as.matrix(results_meta[, prob_cols, drop = FALSE])
      prob_matrix[is.na(prob_matrix)] <- 0
      confidence_scores <- apply(prob_matrix, 1, max, na.rm = TRUE)
    } else {
      confidence_scores <- rep(0, ncol(sce_test))
    }
    
    # Handle NAs
    preds[is.na(preds)] <- "Unknown"
    confidence_scores[!is.finite(confidence_scores)] <- 0
    
    return(list(
      predictions = as.character(preds),
      true_labels = as.character(colData(sce_test)$Ground_Truth_Celltype),
      confidence_scores = as.numeric(confidence_scores),
      cell_ids = colnames(sce_test),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    ))
    
  }, error = function(e) {
    warning(paste("scAnnotatR error:", e$message))
    return(default_return())
  })
}

# Backward compatibility alias
run_scAnnotatR <- run_scAnnotatR_function