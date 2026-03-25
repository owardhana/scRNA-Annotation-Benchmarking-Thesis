







run_CaSTLe_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required libraries
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package not available. Please install from Bioconductor.")
  }
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost package not available. Please install it first: install.packages('xgboost')")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package not available. Please install it first: install.packages('igraph')")
  }
  
  library(SingleCellExperiment)
  library(xgboost)
  library(Seurat)
  
  # Fixed parameters
  BREAKS <- c(-1, 0, 1, 6, Inf)
  nFeatures <- 100
  
  # Default return function for error handling
  default_return <- function() {
    n_test_cells <- ncol(seurat_test)
    return(list(
      predictions = as.character(rep("Unknown", n_test_cells)),
      true_labels = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores = as.numeric(rep(0, n_test_cells)),
      cell_ids = as.character(colnames(seurat_test)),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }

  # Validate input data
  if (!("Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data))) {
    warning("Ground_Truth_Celltype not found in training data")
    return(default_return())
  }
  if (!("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data))) {
    warning("Ground_Truth_Celltype not found in test data")
    return(default_return())
  }
  
  tryCatch({
    cat("Extracting raw counts for CaSTLe...\n")
    
    # Extract raw count data
    counts_train <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
    counts_test <- GetAssayData(seurat_test, assay = "RNA", layer = "counts")
    
    if (is.null(counts_train) || is.null(counts_test)) {
      warning("Failed to extract count data from Seurat objects")
      return(default_return())
    }
    
    cat(sprintf("Training data: %d genes x %d cells\n", nrow(counts_train), ncol(counts_train)))
    
    # Extract cell type labels
    sourceCellTypes <- as.factor(seurat_train$Ground_Truth_Celltype)
    if (any(is.na(sourceCellTypes))) {
      warning("NA values found in training labels")
      return(default_return())
    }
    
    # Step 1: Find common genes
    cat("Finding common genes...\n")
    source_n_cells_counts <- apply(counts_train, 1, function(x) { sum(x > 0) })
    target_n_cells_counts <- apply(counts_test, 1, function(x) { sum(x > 0) })
    
    common_genes <- intersect(
      rownames(counts_train)[source_n_cells_counts > 10],
      rownames(counts_test)[target_n_cells_counts > 10]
    )
    
    if (length(common_genes) < 50) {
      warning(sprintf("Too few common genes found: %d (need at least 50)", length(common_genes)))
      return(default_return())
    }
    
    # Combine datasets
    ds1 <- t(as.matrix(counts_train[common_genes, ])) 
    ds2 <- t(as.matrix(counts_test[common_genes, ]))   
    ds <- rbind(ds1, ds2)
    isSource <- c(rep(TRUE, nrow(ds1)), rep(FALSE, nrow(ds2)))
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Function to run the core pipeline (reused for benchmarking)
    run_castle_core <- function() {
      # Feature selection
      topFeaturesAvg <- colnames(ds)[order(apply(ds, 2, mean), decreasing = TRUE)]
      
      mi_scores <- apply(ds[isSource, ], 2, function(x) {
        tryCatch({
          binned <- cut(x, breaks = BREAKS, include.lowest = TRUE)
          igraph::compare(binned, sourceCellTypes, method = "nmi")
        }, error = function(e) 0)
      })
      topFeaturesMi <- names(sort(mi_scores, decreasing = TRUE))
      
      selectedFeatures <- union(head(topFeaturesAvg, nFeatures), head(topFeaturesMi, nFeatures))
      
      # Correlation filtering
      cor_matrix <- cor(ds[, selectedFeatures], method = "pearson")
      cor_matrix[!lower.tri(cor_matrix)] <- 0
      selectedFeatures <- selectedFeatures[apply(cor_matrix, 2, function(x) any(x < 0.9))]
      
      if (length(selectedFeatures) < 10) stop("Too few features after filtering")
      
      # Binning and Dummy Variables
      dsBins <- apply(ds[, selectedFeatures], 2, cut, breaks = BREAKS, include.lowest = TRUE)
      nUniq <- apply(dsBins, 2, function(x) length(unique(x)))
      dsBins <- dsBins[, nUniq > 1, drop = FALSE]
      
      if (ncol(dsBins) < 5) stop("Too few variable features after binning")
      
      ds_processed <- model.matrix(~ ., as.data.frame(dsBins))
      ds_processed <- ds_processed[, -1, drop = FALSE]
      
      train_data <- ds_processed[isSource, , drop = FALSE]
      test_data <- ds_processed[!isSource, , drop = FALSE]
      
      cell_type_levels <- levels(sourceCellTypes)
      train_labels <- as.numeric(sourceCellTypes) - 1
      
      # --- UPDATED XGBOOST TRAINING ---
      # Create DMatrix objects (required for proper param handling in new XGBoost)
      dtrain <- xgb.DMatrix(data = train_data, label = train_labels)
      dtest <- xgb.DMatrix(data = test_data)
      
      params <- list(
        max_depth = 5,
        learning_rate = 0.7,      # Renamed from 'eta'
        min_split_loss = 0.001,   # Renamed from 'gamma'
        min_child_weight = 10,
        nthread = 5
      )
      
      if (length(cell_type_levels) > 2) {
        # Multiclass
        params$objective <- "multi:softprob"
        params$num_class <- length(cell_type_levels)
      } else {
        # Binary
        params$objective <- "binary:logistic"
      }
      
      # Use xgb.train instead of xgboost() wrapper
      xg_model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = 20,
        verbose = 0
      )
      
      # Prediction
      if (length(cell_type_levels) > 2) {
        pred_probs <- predict(xg_model, dtest, reshape = TRUE)
        pred_class_numeric <- apply(pred_probs, 1, which.max) - 1
        confidence_scores <- apply(pred_probs, 1, max)
      } else {
        pred_probs <- predict(xg_model, dtest)
        pred_class_numeric <- ifelse(pred_probs > 0.5, 1, 0)
        confidence_scores <- ifelse(pred_class_numeric == 1, pred_probs, 1 - pred_probs)
      }
      
      return(list(
        pred_class_numeric = pred_class_numeric,
        confidence_scores = confidence_scores,
        cell_type_levels = cell_type_levels
      ))
    }
    
    # Run with memory tracking if peakRAM is available, otherwise run directly
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      results <- run_castle_core()
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({ results <- run_castle_core() })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }
    
    # Unpack results
    pred_class_numeric <- results$pred_class_numeric
    confidence_scores <- results$confidence_scores
    cell_type_levels <- results$cell_type_levels
    
    # Map predictions
    predicted_labels <- cell_type_levels[pred_class_numeric + 1]
    true_labels <- seurat_test$Ground_Truth_Celltype
    
    # Handle Unknowns
    cell_predictions <- as.character(predicted_labels)
    cell_predictions[is.na(cell_predictions)] <- "Unknown"
    confidence_scores[is.na(confidence_scores)] <- 0
    
    # Summary
    cat("\n=== CaSTLe Prediction Summary ===\n")
    cat(sprintf("Accuracy: %.2f%%\n", 100 * mean(cell_predictions == true_labels, na.rm = TRUE)))
    cat("=================================\n")
    
    return(list(
      predictions = as.character(cell_predictions),
      true_labels = as.character(true_labels),
      confidence_scores = as.numeric(confidence_scores),
      cell_ids = as.character(colnames(seurat_test)),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    ))
    
  }, error = function(e) {
    warning("CaSTLe error: ", e$message)
    return(default_return())
  })
}

run_CaSTLe <- run_CaSTLe_function