# run_scID.R
#################################################
# scID Function for Benchmarking Framework
# Minimal implementation
#################################################

run_scID_function <- function(seurat_train, seurat_test, markers) {

  library(scID)
  library(Seurat)
  library(MAST)

  # Default return for errors
  default_return <- function() {
    list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    )
  }

  tryCatch({

    # Extract RAW COUNTS from Seurat objects (with v4/v5 compatibility)
    train_data <- tryCatch({
      seurat_train@assays$RNA@layers$counts  # Seurat v5
    }, error = function(e) {
      seurat_train@assays$RNA@counts  # Seurat v4
    })

    test_data <- tryCatch({
      seurat_test@assays$RNA@layers$counts  # Seurat v5
    }, error = function(e) {
      seurat_test@assays$RNA@counts  # Seurat v4
    })

    # Set rownames and colnames
    rownames(train_data) <- rownames(seurat_train)
    colnames(train_data) <- colnames(seurat_train)
    rownames(test_data) <- rownames(seurat_test)
    colnames(test_data) <- colnames(seurat_test)

    # Convert to dense matrix then data frame
    # scID expects data frames with genes (rows) x cells (columns) format
    train_mat <- as.matrix(train_data)
    test_mat <- as.matrix(test_data)

    # Create data frames (scID requires data.frame input)
    reference_gem <- data.frame(train_mat, check.names = FALSE, stringsAsFactors = FALSE)
    target_gem <- data.frame(test_mat, check.names = FALSE, stringsAsFactors = FALSE)

    # Create reference_clusters as plain character vector with names
    ref_labels <- as.character(seurat_train$Ground_Truth_Celltype)
    ref_cells <- as.character(colnames(seurat_train))
    reference_clusters <- ref_labels
    names(reference_clusters) <- ref_cells

    # Verify dimensions match
    cat(sprintf("Reference: %d genes x %d cells\n", nrow(reference_gem), ncol(reference_gem)))
    cat(sprintf("Target: %d genes x %d cells\n", nrow(target_gem), ncol(target_gem)))
    cat(sprintf("Reference clusters: %d cells\n", length(reference_clusters)))

    if (ncol(reference_gem) != length(reference_clusters)) {
      stop(sprintf("Mismatch: reference_gem has %d cells but reference_clusters has %d",
                   ncol(reference_gem), length(reference_clusters)))
    }

    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Run scID with memory tracking
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      scID_output <- scid_multiclass(
        target_gem = target_gem,
        reference_gem = reference_gem,
        reference_clusters = reference_clusters,
        logFC = 0.5, #default = 0.5
        normalize_reference = TRUE,  # Let scID normalize raw counts
        estimate_weights_from_target = FALSE,
        only_pos = FALSE
      )
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        scID_output <- scid_multiclass(
          target_gem = target_gem,
          reference_gem = reference_gem,
          reference_clusters = reference_clusters,
          logFC = 0.5, #default = 0.5
          normalize_reference = TRUE,  # Let scID normalize raw counts
          estimate_weights_from_target = FALSE,
          only_pos = FALSE
        )
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    # Extract results
    predictions <- scID_output$labels
    predictions[is.na(predictions) | predictions == "unassigned"] <- "Unknown"

    # Get confidence from scores matrix
    if (!is.null(scID_output$scores)) {
      confidence_scores <- apply(scID_output$scores, 2, max, na.rm = TRUE)
      confidence_scores[is.infinite(confidence_scores)] <- 0
    } else {
      confidence_scores <- rep(1.0, length(predictions))
    }
    confidence_scores[predictions == "Unknown"] <- 0

    # Return
    list(
      predictions = as.character(predictions),
      true_labels = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores = as.numeric(confidence_scores),
      cell_ids = as.character(colnames(seurat_test)),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    )

  }, error = function(e) {
    warning("scID error: ", e$message)
    return(default_return())
  })
}

# Backward compatibility
run_scID <- run_scID_function
