# run_scibetR.R
#################################################
# scibetR Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scibetR Cell Type Annotation Function
#' 
#' Purpose: Run scibetR algorithm using reference-based classification with TPM conversion
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scibetR with TPM conversion and feature selection
run_scibetR_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required libraries
  suppressMessages(library(ggplot2))
  suppressMessages(library(tidyverse))
  suppressMessages(library(scibetR))
  suppressMessages(library(viridis))
  suppressMessages(library(ggsci))
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

  # Run scibetR pipeline with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
    # 1) Extract expression matrix (genes x cells)
    # Use normalized values if you already have TPM; otherwise compute TPM from counts
    counts_ref <- as.matrix(seurat_train@assays$RNA@layers[["data"]])    # genes x cells
    counts_test <- as.matrix(seurat_test@assays$RNA@layers[["data"]])    # genes x cells

    # 2) Convert counts -> TPM per cell (if SciBet expects TPM; user said TPM)
    # counts are genes x cells; we want TPM with same shape
    counts_to_tpm <- function(cnt_mat){
      # avoid dividing by zero:
      libsize <- colSums(cnt_mat)
      libsize[libsize == 0] <- 1
      t( t(cnt_mat) / libsize * 1e6 )
    }

    tpm_ref  <- counts_to_tpm(counts_ref)   # genes x cells
    tpm_test <- counts_to_tpm(counts_test)

    # 3) Transpose to cells x genes (SciBet expects rows = cells)
    expr_ref  <- as.data.frame(t(tpm_ref))   # rows = cells, cols = genes
    expr_test <- as.data.frame(t(tpm_test))

    # 4) Add label column (last column) to reference
    expr_ref$label <- seurat_train$Ground_Truth_Celltype  # ensure same length/order

    # 5) Make sure gene sets (columns) match between train and test
    common_genes <- intersect(colnames(expr_ref)[colnames(expr_ref) != "label"],
                              colnames(expr_test))
    expr_ref2  <- expr_ref[, c(common_genes, "label")]
    expr_test2 <- expr_test[, common_genes]

    # 6) (Optional) Feature selection using SelectGene on the reference
    selected <- SelectGene_R(expr_ref2, k = 1000)
    # Convert selected to gene names; SelectGene may return char vector or logical/indices
    if (is.character(selected)) {
      sel_genes <- selected
    } else if (is.logical(selected)) {
      sel_genes <- names(expr_ref2)[which(selected)]
    } else {
      sel_genes <- names(expr_ref2)[selected]
    }

    # Subset both train/test to the selected genes
    train_set <- expr_ref2[, c(sel_genes, "label")]
    test_set  <- expr_test2[, sel_genes]

    # 7) Run SciBet (train/reference first, test without label column)
    prd <- SciBet_R(train_set, test_set, k = 1000)
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        # 1) Extract expression matrix (genes x cells)
        counts_ref <- as.matrix(seurat_train@assays$RNA@layers[["data"]])
        counts_test <- as.matrix(seurat_test@assays$RNA@layers[["data"]])

        # 2) Convert counts -> TPM per cell
        counts_to_tpm <- function(cnt_mat){
          libsize <- colSums(cnt_mat)
          libsize[libsize == 0] <- 1
          t( t(cnt_mat) / libsize * 1e6 )
        }

        tpm_ref  <- counts_to_tpm(counts_ref)
        tpm_test <- counts_to_tpm(counts_test)

        # 3) Transpose to cells x genes
        expr_ref  <- as.data.frame(t(tpm_ref))
        expr_test <- as.data.frame(t(tpm_test))

        # 4) Add label column to reference
        expr_ref$label <- seurat_train$Ground_Truth_Celltype

        # 5) Match gene sets between train and test
        common_genes <- intersect(colnames(expr_ref)[colnames(expr_ref) != "label"],
                                  colnames(expr_test))
        expr_ref2  <- expr_ref[, c(common_genes, "label")]
        expr_test2 <- expr_test[, common_genes]

        # 6) Feature selection
        selected <- SelectGene_R(expr_ref2, k = 1000)
        if (is.character(selected)) {
          sel_genes <- selected
        } else if (is.logical(selected)) {
          sel_genes <- names(expr_ref2)[which(selected)]
        } else {
          sel_genes <- names(expr_ref2)[selected]
        }

        # Subset both train/test to selected genes
        train_set <- expr_ref2[, c(sel_genes, "label")]
        test_set  <- expr_test2[, sel_genes]

        # 7) Run SciBet
        prd <- SciBet_R(train_set, test_set, k = 1000)
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  # 8) Extract results
  cell_predictions <- prd
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Create confidence scores (SciBet doesn't provide these directly)
  confidence_scores <- rep(1.0, length(cell_predictions))  # Assume high confidence for all predictions

  # Handle NA predictions
  cell_predictions[is.na(cell_predictions)] <- "Unknown"

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
    warning(sprintf("scibetR failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scibetR <- run_scibetR_function