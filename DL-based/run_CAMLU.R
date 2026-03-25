# run_CAMLU.R
#################################################
# CAMLU Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' CAMLU Cell Type Annotation Function
#'
#' Purpose: Run CAMLU algorithm using autoencoder-based method for unknown cell detection
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference for model training)
#'   - seurat_test: Test Seurat object to predict
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses CAMLU's autoencoder with iterative feature selection for novel cell detection
run_CAMLU_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("CAMLU", quietly = TRUE)) {
    stop("CAMLU package not available. Please install it first: devtools::install_github('ziyili20/CAMLU', build_vignettes = FALSE)")
  }
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("keras package not available. Please install it first: install.packages('keras')")
  }

  library(CAMLU)
  library(keras)
  library(Seurat)

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

    cat("Extracting count matrices for CAMLU...\n")

    # Extract count matrices (CAMLU expects raw counts)
    x_train <- tryCatch({
      # Try to get raw counts first
      GetAssayData(seurat_train, assay = "RNA", layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(seurat_train, assay = "RNA", layer = "counts")
      }, error = function(e2) {
        # Fallback to normalized data if raw counts not available
        warning("Raw counts not available, using normalized data")
        GetAssayData(seurat_train, assay = "RNA", layer = "data")
      })
    })

    x_test <- tryCatch({
      # Try to get raw counts first
      GetAssayData(seurat_test, assay = "RNA", layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(seurat_test, assay = "RNA", layer = "counts")
      }, error = function(e2) {
        # Fallback to normalized data if raw counts not available
        GetAssayData(seurat_test, assay = "RNA", layer = "data")
      })
    })

    # Convert sparse matrices to regular matrices if needed
    if (inherits(x_train, "sparseMatrix")) {
      x_train <- as.matrix(x_train)
    }
    if (inherits(x_test, "sparseMatrix")) {
      x_test <- as.matrix(x_test)
    }

    # Get common genes between train and test
    common_genes <- intersect(rownames(x_train), rownames(x_test))

    if (length(common_genes) < 100) {
      warning("Too few common genes between train and test datasets")
      return(default_return())
    }

    # Subset to common genes
    x_train <- x_train[common_genes, ]
    x_test <- x_test[common_genes, ]

    cat(sprintf("Using %d common genes for CAMLU analysis\n", length(common_genes)))

    # Extract training labels for full annotation mode
    y_train <- as.character(seurat_train$Ground_Truth_Celltype)

    # Check if we have sufficient training data
    cell_type_counts <- table(y_train)
    min_cells_per_type <- 3  # Minimum for autoencoder training

    if (any(cell_type_counts < min_cells_per_type)) {
      rare_types <- names(cell_type_counts)[cell_type_counts < min_cells_per_type]
      warning(sprintf("Some cell types have fewer than %d cells: %s",
                      min_cells_per_type, paste(rare_types, collapse = ", ")))
    }

    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells_per_type]
    if (length(valid_cell_types) < 2) {
      warning("Insufficient cell types with enough cells for CAMLU training")
      return(default_return())
    }

    cat(sprintf("Training CAMLU with %d cell types across %d training cells\n",
                length(unique(y_train)), ncol(x_train)))

    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Run CAMLU with full annotation mode and memory tracking
    cat("Running CAMLU autoencoder training and prediction...\n")

    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      # Using full_annotation = TRUE to get complete cell type predictions
      label_result <- CAMLU(
        x_train = x_train,
        x_test = x_test,
        y_train = y_train,
        full_annotation = TRUE,  # Enable full cell type annotation
        ngene = 5000,
        lognormalize = TRUE
      )
      cell_predictions <- label_result$label_full
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        # Using full_annotation = TRUE to get complete cell type predictions
        label_result <- CAMLU(
          x_train = x_train,
          x_test = x_test,
          y_train = y_train,
          full_annotation = TRUE,  # Enable full cell type annotation
          ngene = 5000,
          lognormalize = TRUE
        )
        cell_predictions <- label_result$label_full
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    # Get true labels for test data
    true_labels <- seurat_test$Ground_Truth_Celltype

    # DIAGNOSTIC OUTPUT
    cat("\n=== CAMLU DIAGNOSTIC OUTPUT ===\n")
    cat("Training cell types (y_train unique):\n")
    print(unique(y_train))
    cat("\nCAMLU raw predictions (first 20):\n")
    print(head(cell_predictions, 20))
    cat("\nCAMLU prediction class:\n")
    print(class(cell_predictions))
    cat("\nUnique CAMLU predictions:\n")
    print(unique(cell_predictions))
    cat("\nTrue labels (first 20):\n")
    print(head(true_labels, 20))
    cat("\nUnique true labels:\n")
    print(unique(true_labels))
    cat("\nDirect match count (before any processing):\n")
    print(sum(cell_predictions == true_labels, na.rm = TRUE))
    cat("===============================\n\n")

    # Ensure predictions have correct length
    if (length(cell_predictions) != ncol(seurat_test)) {
      warning(sprintf("Prediction length mismatch. Expected: %d, Got: %d",
                      ncol(seurat_test), length(cell_predictions)))
      return(default_return())
    }

    # Handle missing or NA predictions
    cell_predictions[is.na(cell_predictions)] <- "Unknown"
    cell_predictions[cell_predictions == ""] <- "Unknown"

    cat("\n=== AFTER NA HANDLING ===\n")
    cat("Predictions that became 'Unknown':\n")
    print(sum(cell_predictions == "Unknown"))
    cat("Remaining unique predictions:\n")
    print(unique(cell_predictions))
    cat("========================\n\n")

    # Create confidence scores
    # CAMLU doesn't provide explicit confidence scores, use heuristics
    confidence_scores <- rep(0.8, length(cell_predictions))  # Default confidence
    confidence_scores[cell_predictions == "Unknown"] <- 0

    # If CAMLU detected novel cells, those should have lower confidence
    known_cell_types <- unique(y_train)
    novel_predictions <- !cell_predictions %in% known_cell_types

    cat("\n=== NOVEL CELL DETECTION ===\n")
    cat("Known cell types from training:\n")
    print(known_cell_types)
    cat("\nPredictions NOT in known types (novel):\n")
    print(unique(cell_predictions[novel_predictions]))
    cat(sprintf("Novel prediction count: %d out of %d\n", sum(novel_predictions), length(cell_predictions)))
    cat("===========================\n\n")

    confidence_scores[novel_predictions] <- 0.3  # Lower confidence for novel types

    cat(sprintf("CAMLU completed successfully:\n"))
    cat(sprintf("  - Total predictions: %d\n", length(cell_predictions)))
    cat(sprintf("  - Known cell types: %d\n", sum(!novel_predictions)))
    cat(sprintf("  - Novel/Unknown predictions: %d\n", sum(novel_predictions)))

    # Final accuracy check before returning
    cat(sprintf("  - Final accuracy (raw): %.2f%%\n", 100 * mean(cell_predictions == true_labels, na.rm = TRUE)))

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
    warning("CAMLU error: ", e$message)

    # Check for common keras/tensorflow issues
    if (grepl("tensorflow|keras", e$message, ignore.case = TRUE)) {
      warning("This appears to be a TensorFlow/Keras error. Please ensure TensorFlow is properly installed.")
      warning("You may need to run: keras::install_tensorflow()")
    }

    return(default_return())
  })
}

# For backward compatibility
run_CAMLU <- run_CAMLU_function