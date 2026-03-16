#' # run_CaSTLe.R
#' #################################################
#' # CaSTLe Function for Benchmarking Framework
#' # Input: Train/test Seurat objects and markers from FindAllMarkers
#' # Output: Standardized results format for CV framework
#' # Note: CaSTLe = Cell Annotation using Supervised Transfer Learning with Expression binning
#' # Reference: Uses XGBoost with feature selection based on mean expression and mutual information
#' #################################################
#' 
#' #' CaSTLe Cell Type Annotation Function
#' #'
#' #' Purpose: Run CaSTLe algorithm using XGBoost with intelligent feature selection
#' #' Inputs:
#' #'   - seurat_train: Training Seurat object (reference data)
#' #'   - seurat_test: Test Seurat object (query data to predict)
#' #'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' #' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' #' Algorithm:
#' #'   - Feature selection: Top genes by mean expression + mutual information
#' #'   - Data transformation: Expression binning + dummy variables
#' #'   - Classification: XGBoost gradient boosting
#' #' Dependencies: SingleCellExperiment, xgboost, igraph, Seurat
#' #' Fixed Parameters:
#' #'   - nFeatures = 100 (top genes to select)
#' #'   - BREAKS = c(-1, 0, 1, 6, Inf) (expression binning thresholds)
#' #'   - XGBoost: eta=0.7, nround=20, gamma=0.001, max_depth=5, min_child_weight=10
#' run_CaSTLe_function <- function(seurat_train, seurat_test, markers) {
#' 
#'   # Load required libraries
#'   if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     stop("SingleCellExperiment package not available. Please install from Bioconductor.")
#'   }
#'   if (!requireNamespace("xgboost", quietly = TRUE)) {
#'     stop("xgboost package not available. Please install it first: install.packages('xgboost')")
#'   }
#'   if (!requireNamespace("igraph", quietly = TRUE)) {
#'     stop("igraph package not available. Please install it first: install.packages('igraph')")
#'   }
#' 
#'   library(SingleCellExperiment)
#'   library(xgboost)
#'   # Note: igraph is NOT loaded with library() because Seurat imports it
#'   # Using igraph::compare() with namespace notation instead to avoid conflicts
#'   library(Seurat)
#' 
#'   # Fixed parameters (following original CaSTLe implementation)
#'   BREAKS <- c(-1, 0, 1, 6, Inf)
#'   nFeatures <- 100
#' 
#'   # Default return function for error handling
#'   default_return <- function() {
#'     n_test_cells <- ncol(seurat_test)
#'     return(list(
#'       predictions = as.character(rep("Unknown", n_test_cells)),
#'       true_labels = as.character(seurat_test$Ground_Truth_Celltype),
#'       confidence_scores = as.numeric(rep(0, n_test_cells)),
#'       cell_ids = as.character(colnames(seurat_test)),
#'       peak_memory_mb = NA
#'     ))
#'   }
#' 
#'   # Validate input data
#'   if (!("Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data))) {
#'     warning("Ground_Truth_Celltype not found in training data")
#'     return(default_return())
#'   }
#' 
#'   if (!("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data))) {
#'     warning("Ground_Truth_Celltype not found in test data")
#'     return(default_return())
#'   }
#' 
#'   tryCatch({
#' 
#'     cat("Extracting raw counts for CaSTLe...\\n")
#' 
#'     # Extract raw count data (CaSTLe uses raw counts, not normalized)
#'     # Use cross-compatible approach for Seurat v4/v5
#'     counts_train <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
#' 
#'     counts_test <- GetAssayData(seurat_test, assay = "RNA", layer = "counts")
#' 
#'     # Validate data extraction
#'     if (is.null(counts_train) || is.null(counts_test)) {
#'       warning("Failed to extract count data from Seurat objects")
#'       return(default_return())
#'     }
#' 
#'     cat(sprintf("Training data: %d genes x %d cells\\n", nrow(counts_train), ncol(counts_train)))
#'     cat(sprintf("Test data: %d genes x %d cells\\n", nrow(counts_test), ncol(counts_test)))
#' 
#'     # Extract cell type labels
#'     sourceCellTypes <- as.factor(seurat_train$Ground_Truth_Celltype)
#' 
#'     # Validate labels
#'     if (any(is.na(sourceCellTypes))) {
#'       warning("NA values found in training labels")
#'       return(default_return())
#'     }
#' 
#'     cat(sprintf("Cell types in training data: %d unique\\n", length(unique(sourceCellTypes))))
#' 
#'     # Step 1: Find common genes (expressed in >10 cells in both datasets)
#'     cat("Finding common genes...\\n")
#'     source_n_cells_counts <- apply(counts_train, 1, function(x) { sum(x > 0) })
#'     target_n_cells_counts <- apply(counts_test, 1, function(x) { sum(x > 0) })
#' 
#'     common_genes <- intersect(
#'       rownames(counts_train)[source_n_cells_counts > 10],
#'       rownames(counts_test)[target_n_cells_counts > 10]
#'     )
#' 
#'     if (length(common_genes) < 50) {
#'       warning(sprintf("Too few common genes found: %d (need at least 50)", length(common_genes)))
#'       return(default_return())
#'     }
#' 
#'     cat(sprintf("Common genes found: %d\\n", length(common_genes)))
#' 
#'     # Subset to common genes
#'     ds1 <- t(as.matrix(counts_train[common_genes, ]))  # cells x genes
#'     ds2 <- t(as.matrix(counts_test[common_genes, ]))   # cells x genes
#' 
#'     # Combine datasets for unified processing
#'     ds <- rbind(ds1, ds2)
#'     isSource <- c(rep(TRUE, nrow(ds1)), rep(FALSE, nrow(ds2)))
#' 
#'     cat(sprintf("Combined dataset: %d cells x %d genes\\n", nrow(ds), ncol(ds)))
#' 
#'     # Track peak memory usage
#'     peak_memory_mb <- NA
#' 
#'     # Run CaSTLe feature selection + training + prediction with memory tracking
#'     if (!requireNamespace("bench", quietly = TRUE)) {
#'       warning("bench package not available for memory tracking")
#' 
#'       # Step 2: Feature selection - Top genes by mean expression
#'       cat("Selecting features by mean expression...\\n")
#'       topFeaturesAvg <- colnames(ds)[order(apply(ds, 2, mean), decreasing = TRUE)]
#' 
#'       # Step 3: Feature selection - Top genes by mutual information with cell types
#'       cat("Calculating mutual information with cell types...\\n")
#' 
#'       # Calculate NMI for each gene in source data
#'       mi_scores <- apply(ds[isSource, ], 2, function(x) {
#'         tryCatch({
#'           # Bin the expression values for MI calculation
#'           binned <- cut(x, breaks = BREAKS, include.lowest = TRUE)
#'           # Calculate normalized mutual information
#'           igraph::compare(binned, sourceCellTypes, method = "nmi")
#'         }, error = function(e) {
#'           return(0)  # Return 0 if calculation fails
#'         })
#'       })
#' 
#'       topFeaturesMi <- names(sort(mi_scores, decreasing = TRUE))
#' 
#'       # Step 4: Union of top n genes from both methods
#'       cat(sprintf("Selecting union of top %d genes from both methods...\\n", nFeatures))
#'       selectedFeatures <- union(head(topFeaturesAvg, nFeatures), head(topFeaturesMi, nFeatures))
#' 
#'       cat(sprintf("Selected features before correlation filtering: %d\\n", length(selectedFeatures)))
#' 
#'       # Step 5: Remove highly correlated features (r > 0.9)
#'       cat("Removing highly correlated features...\\n")
#'       cor_matrix <- cor(ds[, selectedFeatures], method = "pearson")
#'       cor_matrix[!lower.tri(cor_matrix)] <- 0
#' 
#'       # Keep features that don't have high correlation with any other feature
#'       selectedFeatures <- selectedFeatures[apply(cor_matrix, 2, function(x) any(x < 0.9))]
#' 
#'       cat(sprintf("Selected features after correlation filtering: %d\\n", length(selectedFeatures)))
#' 
#'       if (length(selectedFeatures) < 10) {
#'         warning("Too few features after filtering")
#'         return(default_return())
#'       }
#' 
#'       # Step 6: Bin expression values into discrete categories
#'       cat("Binning expression values...\\n")
#'       dsBins <- apply(ds[, selectedFeatures], 2, cut, breaks = BREAKS, include.lowest = TRUE)
#' 
#'       # Step 7: Use only bins with more than one unique value
#'       nUniq <- apply(dsBins, 2, function(x) { length(unique(x)) })
#'       dsBins <- dsBins[, nUniq > 1, drop = FALSE]
#' 
#'       cat(sprintf("Features with sufficient variation: %d\\n", ncol(dsBins)))
#' 
#'       if (ncol(dsBins) < 5) {
#'         warning("Too few variable features after binning")
#'         return(default_return())
#'       }
#' 
#'       # Step 8: Convert bins to dummy variables
#'       cat("Converting to dummy variables...\\n")
#'       ds_processed <- model.matrix(~ ., as.data.frame(dsBins))
#' 
#'       # Remove intercept column
#'       ds_processed <- ds_processed[, -1, drop = FALSE]
#' 
#'       cat(sprintf("Final feature matrix: %d cells x %d features\\n", nrow(ds_processed), ncol(ds_processed)))
#' 
#'       # Step 9: Prepare training data (100% of source data, no internal validation split)
#'       cat("Training XGBoost classifier...\\n")
#'       train_data <- ds_processed[isSource, , drop = FALSE]
#'       test_data <- ds_processed[!isSource, , drop = FALSE]
#' 
#'       # Convert cell type labels to numeric (0-indexed for xgboost)
#'       cell_type_levels <- levels(sourceCellTypes)
#'       train_labels <- as.numeric(sourceCellTypes) - 1
#' 
#'       # Step 10: Train XGBoost
#'       # Use different objective for multiclass vs binary classification
#'       if (length(cell_type_levels) > 2) {
#'         # Multiclass classification
#'         xg_model <- xgboost(
#'           data = train_data,
#'           label = train_labels,
#'           objective = "multi:softprob",  # Use softprob to get probabilities
#'           num_class = length(cell_type_levels),
#'           eta = 0.7,
#'           nthread = 5,
#'           nround = 20,
#'           verbose = 0,
#'           gamma = 0.001,
#'           max_depth = 5,
#'           min_child_weight = 10
#'         )
#'       } else {
#'         # Binary classification
#'         xg_model <- xgboost(
#'           data = train_data,
#'           label = train_labels,
#'           objective = "binary:logistic",
#'           eta = 0.7,
#'           nthread = 5,
#'           nround = 20,
#'           verbose = 0,
#'           gamma = 0.001,
#'           max_depth = 5,
#'           min_child_weight = 10
#'         )
#'       }
#' 
#'       cat("XGBoost training completed.\\n")
#' 
#'       # Step 11: Predict on test data
#'       cat("Predicting cell types for test data...\\n")
#' 
#'       if (length(cell_type_levels) > 2) {
#'         # Multiclass: predict returns probability matrix
#'         pred_probs <- predict(xg_model, test_data, reshape = TRUE)
#' 
#'         # Get predicted class (highest probability)
#'         pred_class_numeric <- apply(pred_probs, 1, which.max) - 1
#' 
#'         # Get confidence scores (max probability)
#'         confidence_scores <- apply(pred_probs, 1, max)
#' 
#'       } else {
#'         # Binary: predict returns probabilities for class 1
#'         pred_probs <- predict(xg_model, test_data)
#' 
#'         # Convert to class predictions (threshold at 0.5)
#'         pred_class_numeric <- ifelse(pred_probs > 0.5, 1, 0)
#' 
#'         # Confidence is the probability of predicted class
#'         confidence_scores <- ifelse(pred_class_numeric == 1, pred_probs, 1 - pred_probs)
#'       }
#'     } else {
#'       library(bench)
#'       bench_result <- bench::mark(
#'         {
#'           # Complete CaSTLe pipeline (steps 2-11)
#'           topFeaturesAvg <- colnames(ds)[order(apply(ds, 2, mean), decreasing = TRUE)]
#' 
#'           mi_scores <- apply(ds[isSource, ], 2, function(x) {
#'             tryCatch({
#'               binned <- cut(x, breaks = BREAKS, include.lowest = TRUE)
#'               igraph::compare(binned, sourceCellTypes, method = "nmi")
#'             }, error = function(e) { return(0) })
#'           })
#' 
#'           topFeaturesMi <- names(sort(mi_scores, decreasing = TRUE))
#'           selectedFeatures <- union(head(topFeaturesAvg, nFeatures), head(topFeaturesMi, nFeatures))
#' 
#'           cor_matrix <- cor(ds[, selectedFeatures], method = "pearson")
#'           cor_matrix[!lower.tri(cor_matrix)] <- 0
#'           selectedFeatures <- selectedFeatures[apply(cor_matrix, 2, function(x) any(x < 0.9))]
#' 
#'           if (length(selectedFeatures) < 10) {
#'             stop("Too few features after filtering")
#'           }
#' 
#'           dsBins <- apply(ds[, selectedFeatures], 2, cut, breaks = BREAKS, include.lowest = TRUE)
#'           nUniq <- apply(dsBins, 2, function(x) { length(unique(x)) })
#'           dsBins <- dsBins[, nUniq > 1, drop = FALSE]
#' 
#'           if (ncol(dsBins) < 5) {
#'             stop("Too few variable features after binning")
#'           }
#' 
#'           ds_processed <- model.matrix(~ ., as.data.frame(dsBins))
#'           ds_processed <- ds_processed[, -1, drop = FALSE]
#' 
#'           train_data <- ds_processed[isSource, , drop = FALSE]
#'           test_data <- ds_processed[!isSource, , drop = FALSE]
#' 
#'           cell_type_levels <- levels(sourceCellTypes)
#'           train_labels <- as.numeric(sourceCellTypes) - 1
#' 
#'           if (length(cell_type_levels) > 2) {
#'             xg_model <- xgboost(
#'               data = train_data,
#'               label = train_labels,
#'               objective = "multi:softprob",
#'               num_class = length(cell_type_levels),
#'               eta = 0.7,
#'               nthread = 5,
#'               nround = 20,
#'               verbose = 0,
#'               gamma = 0.001,
#'               max_depth = 5,
#'               min_child_weight = 10
#'             )
#'           } else {
#'             xg_model <- xgboost(
#'               data = train_data,
#'               label = train_labels,
#'               objective = "binary:logistic",
#'               eta = 0.7,
#'               nthread = 5,
#'               nround = 20,
#'               verbose = 0,
#'               gamma = 0.001,
#'               max_depth = 5,
#'               min_child_weight = 10
#'             )
#'           }
#' 
#'           if (length(cell_type_levels) > 2) {
#'             pred_probs <- predict(xg_model, test_data, reshape = TRUE)
#'             pred_class_numeric <- apply(pred_probs, 1, which.max) - 1
#'             confidence_scores <- apply(pred_probs, 1, max)
#'           } else {
#'             pred_probs <- predict(xg_model, test_data)
#'             pred_class_numeric <- ifelse(pred_probs > 0.5, 1, 0)
#'             confidence_scores <- ifelse(pred_class_numeric == 1, pred_probs, 1 - pred_probs)
#'           }
#'         },
#'         memory = TRUE,
#'         iterations = 1,
#'         check = FALSE
#'       )
#'       peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
#' 
#'       # Check for early termination
#'       if (!exists("pred_class_numeric") || !exists("confidence_scores") || !exists("cell_type_levels")) {
#'         warning("CaSTLe pipeline failed")
#'         return_val <- default_return()
#'         return_val$peak_memory_mb <- peak_memory_mb
#'         return(return_val)
#'       }
#'     }
#' 
#'     # Map numeric predictions back to cell type names
#'     predicted_labels <- cell_type_levels[pred_class_numeric + 1]
#' 
#'     # Get true labels for test data
#'     true_labels <- seurat_test$Ground_Truth_Celltype
#' 
#'     # Validate prediction length
#'     if (length(predicted_labels) != ncol(seurat_test)) {
#'       warning(sprintf("Prediction length mismatch. Expected: %d, Got: %d",
#'                       ncol(seurat_test), length(predicted_labels)))
#'       return(default_return())
#'     }
#' 
#'     # Handle missing or NA predictions
#'     cell_predictions <- as.character(predicted_labels)
#'     cell_predictions[is.na(cell_predictions)] <- "Unknown"
#'     cell_predictions[cell_predictions == ""] <- "Unknown"
#'     confidence_scores[is.na(confidence_scores)] <- 0
#' 
#'     # Print prediction summary
#'     cat("\\n=== CaSTLe Prediction Summary ===\\n")
#'     cat("Unique predictions:\\n")
#'     print(table(cell_predictions))
#'     cat(sprintf("Total predictions: %d\\n", length(cell_predictions)))
#'     cat(sprintf("Unknown predictions: %d\\n", sum(cell_predictions == "Unknown")))
#'     cat(sprintf("Mean confidence: %.3f\\n", mean(confidence_scores, na.rm = TRUE)))
#'     cat(sprintf("Accuracy: %.2f%%\\n", 100 * mean(cell_predictions == true_labels, na.rm = TRUE)))
#'     cat("=================================\\n\\n")
#' 
#'     # Return standardized format
#'     return(list(
#'       predictions = as.character(cell_predictions),
#'       true_labels = as.character(true_labels),
#'       confidence_scores = as.numeric(confidence_scores),
#'       cell_ids = as.character(colnames(seurat_test)),
#'       peak_memory_mb = peak_memory_mb
#'     ))
#' 
#'   }, error = function(e) {
#'     warning("CaSTLe error: ", e$message)
#' 
#'     # Check for common package issues
#'     if (grepl("igraph|compare", e$message, ignore.case = TRUE)) {
#'       warning("This appears to be an igraph error. Please ensure igraph is properly installed.")
#'       warning("You may need to run: install.packages('igraph')")
#'     }
#' 
#'     if (grepl("xgboost", e$message, ignore.case = TRUE)) {
#'       warning("This appears to be an xgboost error. Please ensure xgboost is properly installed.")
#'       warning("You may need to run: install.packages('xgboost')")
#'     }
#' 
#'     return(default_return())
#'   })
#' }
#' 
#' # For backward compatibility
#' run_CaSTLe <- run_CaSTLe_function











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
      peak_memory_mb = NA
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
    peak_memory_mb <- NA
    
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
    
    # Run with memory tracking if bench is available, otherwise run directly
    if (requireNamespace("bench", quietly = TRUE)) {
      bench_result <- bench::mark(
        results <- run_castle_core(),
        memory = TRUE, iterations = 1, check = FALSE
      )
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
    } else {
      results <- run_castle_core()
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
      peak_memory_mb = peak_memory_mb
    ))
    
  }, error = function(e) {
    warning("CaSTLe error: ", e$message)
    return(default_return())
  })
}

run_CaSTLe <- run_CaSTLe_function