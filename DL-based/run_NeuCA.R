# run_NeuCA.R
#################################################
# NeuCA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
# Note: NeuCA package is deprecated - using local source files from DL-based/NeuCA_helpers/
# Implementation: Bypasses SingleCellExperiment to avoid BiocGenerics corruption issues
#################################################

#' NeuCA Cell Type Annotation Function
#'
#' Purpose: Run NeuCA algorithm using neural networks with hierarchical cell type detection
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference for model training, normalized/log-transformed)
#'   - seurat_test: Test Seurat object to predict (normalized/log-transformed)
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency, not directly used)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm:
#'   - Uses NeuCA's neural network approach with automatic route selection
#'   - Route 1 (Direct NN): If all cell type correlations < 0.95
#'   - Route 2 (Hierarchical NN): If high correlation detected (>= 0.95)
#'   - Model size: "medium" (128 -> 64 -> nClass architecture)
#' Dependencies: keras, e1071 (NO limma/statmod/SingleCellExperiment/BiocGenerics required)
#' Source Files: DL-based/NeuCA_helpers/{functions_NN.R, functions_mhNN.R, utils.R}
#' Implementation: Calls NeuCA internal functions directly with matrices (bypasses SCE for compatibility)
#' Note: limma replaced with simple t-test in GetFeature() to avoid statmod corruption
run_NeuCA_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  # Note: SingleCellExperiment and limma NOT required - using direct NeuCA internal calls
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("keras package not available. Please install it first: install.packages('keras')")
  }
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("e1071 package not available. Please install it first: install.packages('e1071')")
  }

  library(keras)
  library(e1071)
  library(Seurat)

  # Source NeuCA helper functions (local source - deprecated package)
  helper_path <- "DL-based/NeuCA_helpers"
  tryCatch({
    source(file.path(helper_path, "functions_NN.R"))
    source(file.path(helper_path, "functions_mhNN.R"))
    source(file.path(helper_path, "utils.R"))
  }, error = function(e) {
    stop("Failed to load NeuCA source files from ", helper_path, ": ", e$message)
  })

  # Verify NeuCA function loaded
  if (!exists("NeuCA")) {
    stop("NeuCA function not found after sourcing helper files. Check DL-based/NeuCA_helpers/ directory")
  }

  #' Direct NeuCA Wrapper - Bypasses SingleCellExperiment
  #'
  #' Calls NeuCA internal functions directly to avoid BiocGenerics corruption
  #' Implements exact same logic as NeuCA() but with matrices instead of SCE
  #'
  #' @param train_matrix Normalized gene expression matrix (genes x cells)
  #' @param test_matrix Normalized gene expression matrix (genes x cells)
  #' @param train_labels Character vector of cell type labels
  #' @param model.size "small", "medium", or "big"
  #' @param verbose Logical, show training progress
  #' @return Character vector of predicted cell types
  run_NeuCA_direct <- function(train_matrix, test_matrix, train_labels,
                                model.size = "medium", verbose = FALSE) {

    message("Working on scRNA-seq data cell label training and testing (direct mode):")

    # Validation
    if (!model.size %in% c("big", "medium", "small")) {
      stop("model.size must be one of: 'big', 'medium', 'small'")
    }

    if (ncol(train_matrix) != length(train_labels)) {
      stop("Number of cells in train_matrix must match length of train_labels")
    }

    # Step 1: Data normalization (same as NeuCA utils.R:41-43)
    message("Normalizing data...")
    normedData <- SampleNorm(
      train_count = train_matrix,  # genes x cells
      test_count = test_matrix,     # genes x cells
      nfeature = 5000
    )

    # Step 2: Determine correlation (same as NeuCA utils.R:46-49)
    message("Checking cell type correlations...")
    train_labels_char <- as.character(train_labels)
    cd <- cor.det(
      dat = t(normedData$train_count_out),  # cells x genes
      lb = train_labels_char
    )

    # Step 3: Create label matrix for neural networks (same as NeuCA utils.R:51-60)
    CTname <- unique(train_labels_char)
    tr.lab <- matrix(0,
                     nrow = nrow(normedData$train_count_out),
                     ncol = length(CTname))

    for(i in seq_along(CTname)) {
      idx <- which(train_labels_char == CTname[i])
      tr.lab[idx, i] <- 1
    }

    # Step 4: Route selection based on correlation
    if (!isTRUE(cd == 1) && !isTRUE(cd == 2)) {
      message("Warning: cor.det returned unexpected value (", cd, "), defaulting to direct NN (route 1)")
      cd <- 1
    }

    if (cd == 1) {
      # Route 1: Direct Neural Network (low correlation)
      message("Based on correlation values, direct neural network IS adopted")

      if (model.size == "big") {
        message("Neural network: big")
        outmodel <- BigNN(
          train_count = normedData$train_count_out,
          train_label = tr.lab,
          lossname = 'binary_crossentropy',
          last_act = "sigmoid",
          nClass = length(CTname),
          verbose = verbose
        )
      } else if (model.size == "medium") {
        message("Neural network: medium")
        outmodel <- MediumNN(
          train_count = normedData$train_count_out,
          train_label = tr.lab,
          lossname = 'binary_crossentropy',
          last_act = "sigmoid",
          nClass = length(CTname),
          verbose = verbose
        )
      } else if (model.size == "small") {
        message("Neural network: small")
        outmodel <- SmallNN(
          train_count = normedData$train_count_out,
          train_label = tr.lab,
          lossname = 'binary_crossentropy',
          last_act = "sigmoid",
          nClass = length(CTname),
          verbose = verbose
        )
      }

      # Predict (same as NeuCA utils.R:95-100)
      preres <- predict(outmodel,
                        x = as.matrix(normedData$test_count_out),
                        batch_size = 256,
                        verbose = verbose)

      prenum <- apply(preres, 1, findm)
      predict.label <- CTname[prenum]

    } else if (cd == 2) {
      # Route 2: Marker-guided hierarchical NN (high correlation)
      message("Based on correlation values, marker-guided hierarchical neural network IS adopted")
      message("Marker-guided hierarchical neural network: ", model.size)

      predict.label <- SCHwrapper(
        train_count = t(normedData$train_count_out),
        train_label = train_labels_char,
        test_count = t(as.matrix(normedData$test_count_out)),
        nMark = 500,
        modeltype = model.size,
        CellTypeMark_thres = 3,
        GeneMarkerExpr_thres = 20,
        verbose = verbose
      )
    }

    return(predict.label)
  }

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
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
    warning("Ground_Truth_Celltype not found in training data")
    return(default_return())
  }

  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in test data")
    return(default_return())
  }

  tryCatch({

    cat("Extracting normalized data for NeuCA...\n")

    # Extract normalized, log-transformed data (already present in Seurat objects)
    # Use cross-compatible approach for Seurat v4/v5
    norm_data_train <- tryCatch({
      # Try Seurat v5 layer approach
      GetAssayData(seurat_train, assay = "RNA", layer = "data")
    }, error = function(e) {
      # Fallback to Seurat v4 slot approach
      GetAssayData(seurat_train, assay = "RNA", layer = "data")
    })

    norm_data_test <- tryCatch({
      # Try Seurat v5 layer approach
      GetAssayData(seurat_test, assay = "RNA", layer = "data")
    }, error = function(e) {
      # Fallback to Seurat v4 slot approach
      GetAssayData(seurat_test, assay = "RNA", layer = "data")
    })

    # Validate data extraction
    if (is.null(norm_data_train) || is.null(norm_data_test)) {
      warning("Failed to extract normalized data from Seurat objects")
      return(default_return())
    }

    # Helper function: Safe sparse-to-dense conversion with memory estimation
    convert_to_dense_safe <- function(sparse_mat, mat_name = "matrix") {
      # If already dense, return as-is
      if (!inherits(sparse_mat, "sparseMatrix")) {
        return(as.matrix(sparse_mat))
      }

      # Estimate memory requirement
      n_elements <- as.numeric(nrow(sparse_mat)) * as.numeric(ncol(sparse_mat))
      estimated_size_gb <- (n_elements * 8) / (1024^3)  # 8 bytes per double

      if (estimated_size_gb > 1.0) {
        cat(sprintf("⚠ Large %s detected: %.2f GB memory required for dense conversion\n",
                    mat_name, estimated_size_gb))
        cat("Converting sparse to dense (this may take a moment)...\n")
      }

      # Convert to dense
      dense_mat <- as.matrix(sparse_mat)
      return(dense_mat)
    }

    # Convert sparse matrices to dense matrices (NeuCA requires dense matrices)
    norm_data_train <- convert_to_dense_safe(norm_data_train, "training matrix")
    norm_data_test <- convert_to_dense_safe(norm_data_test, "test matrix")

    cat(sprintf("Training data: %d genes x %d cells\n", nrow(norm_data_train), ncol(norm_data_train)))
    cat(sprintf("Test data: %d genes x %d cells\n", nrow(norm_data_test), ncol(norm_data_test)))

    # Get common genes between train and test
    common_genes <- intersect(rownames(norm_data_train), rownames(norm_data_test))

    if (length(common_genes) < 100) {
      warning("Too few common genes between train and test datasets")
      return(default_return())
    }

    # Subset to common genes
    norm_data_train <- norm_data_train[common_genes, ]
    norm_data_test <- norm_data_test[common_genes, ]

    # Optional: Pre-filter genes for very large datasets to reduce memory
    # NeuCA's SampleNorm() will further filter to top 5000 genes anyway
    if (length(common_genes) > 10000) {
      cat(sprintf("Large gene set detected (%d genes)\n", length(common_genes)))
      cat("Pre-filtering to top 10000 most variable genes before SCE conversion...\n")

      # Calculate gene variance on training data
      gene_vars <- apply(norm_data_train, 1, var)
      top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:10000])

      # Subset both datasets
      norm_data_train <- norm_data_train[top_genes, ]
      norm_data_test <- norm_data_test[top_genes, ]

      cat(sprintf("Reduced to %d genes\n", nrow(norm_data_train)))
    }

    cat(sprintf("Using %d common genes for NeuCA analysis\n", nrow(norm_data_train)))

    # Prepare data for direct NeuCA call (bypassing SingleCellExperiment to avoid BiocGenerics corruption)
    cat("Preparing matrices for direct NeuCA call...\n")

    # Extract cell type labels
    train_labels <- seurat_train$Ground_Truth_Celltype

    # Validate labels
    if (any(is.na(train_labels))) {
      warning("NA values found in training labels")
      return(default_return())
    }

    # Display dataset info
    cat(sprintf("Dataset preparation complete:\n"))
    cat(sprintf("  Training: %d genes x %d cells\n", nrow(norm_data_train), ncol(norm_data_train)))
    cat(sprintf("  Test: %d genes x %d cells\n", nrow(norm_data_test), ncol(norm_data_test)))
    cat(sprintf("  Cell types: %d unique\n", length(unique(train_labels))))

    cell_type_counts <- table(train_labels)
    cat(sprintf("Training cell type distribution:\n"))
    print(cell_type_counts)

    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Run NeuCA using direct internal function calls with memory tracking
    cat("\nRunning NeuCA neural network training and prediction (direct mode)...\n")

    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      predicted_labels <- run_NeuCA_direct(
        train_matrix = norm_data_train,  # genes x cells
        test_matrix = norm_data_test,    # genes x cells
        train_labels = train_labels,
        model.size = "medium",  # Fixed model size for consistency
        verbose = FALSE
      )
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        predicted_labels <- run_NeuCA_direct(
          train_matrix = norm_data_train,  # genes x cells
          test_matrix = norm_data_test,    # genes x cells
          train_labels = train_labels,
          model.size = "medium",  # Fixed model size for consistency
          verbose = FALSE
        )
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    cat("NeuCA prediction completed.\n")

    # Get true labels for test data
    true_labels <- seurat_test$Ground_Truth_Celltype

    # Validate prediction length
    if (length(predicted_labels) != ncol(seurat_test)) {
      warning(sprintf("Prediction length mismatch. Expected: %d, Got: %d",
                      ncol(seurat_test), length(predicted_labels)))
      return(default_return())
    }

    # Handle missing or NA predictions
    cell_predictions <- as.character(predicted_labels)
    cell_predictions[is.na(cell_predictions)] <- "Unknown"
    cell_predictions[cell_predictions == ""] <- "Unknown"

    # Create confidence scores (NeuCA doesn't provide explicit confidence)
    # Use heuristic approach similar to other DL-based tools
    confidence_scores <- rep(0.8, length(cell_predictions))  # Default confidence
    confidence_scores[cell_predictions == "Unknown"] <- 0

    # Print prediction summary
    cat("\n=== NeuCA Prediction Summary ===\n")
    cat("Unique predictions:\n")
    print(table(cell_predictions))
    cat(sprintf("Total predictions: %d\n", length(cell_predictions)))
    cat(sprintf("Unknown predictions: %d\n", sum(cell_predictions == "Unknown")))
    cat(sprintf("Accuracy: %.2f%%\n", 100 * mean(cell_predictions == true_labels, na.rm = TRUE)))
    cat("================================\n\n")

    # Return standardized format
    return(list(
      predictions = as.character(cell_predictions),
      true_labels = as.character(true_labels),
      confidence_scores = as.numeric(confidence_scores),
      cell_ids = as.character(colnames(seurat_test)),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    ))

  }, error = function(e) {
    warning("NeuCA error: ", e$message)

    # Check for common keras/tensorflow issues
    if (grepl("tensorflow|keras", e$message, ignore.case = TRUE)) {
      warning("This appears to be a TensorFlow/Keras error. Please ensure TensorFlow is properly installed.")
      warning("You may need to run: keras::install_tensorflow()")
    }

    return(default_return())
  })
}

# For backward compatibility
run_NeuCA <- run_NeuCA_function
