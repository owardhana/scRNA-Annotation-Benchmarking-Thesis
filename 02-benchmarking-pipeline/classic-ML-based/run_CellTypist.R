# run_CellTypist.R
#################################################
#  CellTypist Implementation using File-based Communication
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#'  CellTypist Cell Type Annotation Function
#' 
#' Purpose: Run CellTypist using file-based communication with Python
#' Inputs:
#'   - seurat_train: Training Seurat object (used to train classifier)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses temporary CSV files and standalone Python script
#' 
#' REQUIREMENTS:
#' - Conda environment 'celltypist_benchmark' with CellTypist installed
#' - Python script 'celltypist_helper.py' in same directory
#' 
#' SETUP:
#' conda create -n celltypist_benchmark python=3.9
#' conda activate celltypist_benchmark
#' conda install -c bioconda -c conda-forge celltypist
run_CellTypist_function <- function(seurat_train, seurat_test, markers) {
  
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
  
  message("🔬 Running CellTypist with file-based approach...")
  
  tryCatch({
    
    # Create temporary directory and unique session ID
    temp_dir <- tempdir()
    session_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1))
    
    message("Session ID: ", session_id)
    message("Temporary directory: ", temp_dir)
    
    # Define file paths with explicit character conversion
    train_file <- as.character(file.path(temp_dir, paste0("train_expr_", session_id, ".csv")))
    test_file <- as.character(file.path(temp_dir, paste0("test_expr_", session_id, ".csv")))
    labels_file <- as.character(file.path(temp_dir, paste0("labels_", session_id, ".txt")))
    genes_file <- as.character(file.path(temp_dir, paste0("genes_", session_id, ".txt")))
    output_file <- as.character(file.path(temp_dir, paste0("predictions_", session_id, ".csv")))
    status_file <- as.character(file.path(temp_dir, paste0("status_", session_id, ".txt")))
    
    # Extract expression matrices (use raw counts)
    message("Extracting expression matrices...")
    
    # Get raw count data with cross-compatible approach
    train_counts <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")

    test_counts <- GetAssayData(seurat_test, assay = "RNA", layer = "counts")
    
    # Convert sparse matrices to regular matrices
    if (inherits(train_counts, "sparseMatrix")) {
      message("Converting training sparse matrix to dense...")
      train_counts <- as.matrix(train_counts)
    }
    
    if (inherits(test_counts, "sparseMatrix")) {
      message("Converting test sparse matrix to dense...")
      test_counts <- as.matrix(test_counts)
    }
    
    # Get common genes between train and test
    common_genes <- intersect(rownames(train_counts), rownames(test_counts))
    
    if (length(common_genes) < 100) {
      warning(paste("Too few common genes:", length(common_genes)))
      return(default_return())
    }
    
    message("Using ", length(common_genes), " common genes")
    
    # Subset to common genes
    train_counts_subset <- train_counts[common_genes, ]
    test_counts_subset <- test_counts[common_genes, ]
    
    # Transpose to cells x genes (CellTypist format)
    train_expr <- t(train_counts_subset)
    test_expr <- t(test_counts_subset)
    
    # Get cell type labels for training data
    train_labels <- seurat_train$Ground_Truth_Celltype
    
    # Validate dimensions
    if (nrow(train_expr) != length(train_labels)) {
      warning("Dimension mismatch between expression matrix and labels")
      return(default_return())
    }
    
    message("Writing temporary files...")
    
    # Write expression matrices (cells x genes)
    write.csv(train_expr, train_file, row.names = TRUE, quote = FALSE)
    write.csv(test_expr, test_file, row.names = TRUE, quote = FALSE)
    
    # Write labels (one per line, corresponding to train cells)
    writeLines(as.character(train_labels), labels_file)
    
    # Write gene names (one per line, corresponding to matrix columns)
    writeLines(common_genes, genes_file)
    
    # Validate files were created
    required_files <- c(train_file, test_file, labels_file, genes_file)
    missing_files <- required_files[!file.exists(required_files)]
    
    if (length(missing_files) > 0) {
      warning("Failed to create temporary files: ", paste(missing_files, collapse = ", "))
      return(default_return())
    }
    
    message("✓ Temporary files created successfully")
    
    # Find Python script path - look in Classic-ML-Based directory
    script_dir <- tryCatch({
      frame_info <- sys.frame(1)
      if ("ofile" %in% names(frame_info) && !is.null(frame_info$ofile)) {
        dirname(as.character(frame_info$ofile))
      } else {
        # Look for Classic-ML-Based directory
        current_dir <- getwd()
        if (grepl("Classic-ML-Based", current_dir)) {
          current_dir
        } else {
          file.path(current_dir, "Classic-ML-Based")
        }
      }
    }, error = function(e) {
      # Default fallback - look for Classic-ML-Based directory
      current_dir <- getwd()
      if (grepl("Classic-ML-Based", current_dir)) {
        current_dir
      } else {
        file.path(current_dir, "Classic-ML-Based")
      }
    })
    
    # Ensure we're looking in the right directory
    if (!grepl("Classic-ML-Based", script_dir)) {
      script_dir <- file.path(script_dir, "Classic-ML-Based")
    }
    
    python_script <- as.character(file.path(script_dir, "celltypist_helper.py"))
    
    if (!file.exists(python_script)) {
      warning("Python script not found: ", python_script)
      return(default_return())
    }
    
    message("Using Python script: ", python_script)
    
    # Find conda environment python path
    # Try common conda locations
    possible_conda_paths <- c(
      "~/anaconda3/envs/celltypist_benchmark/bin/python",
      "~/miniconda3/envs/celltypist_benchmark/bin/python", 
      "~/miniforge3/envs/celltypist_benchmark/bin/python",
      "/opt/anaconda3/envs/celltypist_benchmark/bin/python",
      "/opt/miniconda3/envs/celltypist_benchmark/bin/python"
    )
    
    conda_python <- NULL
    for (path in possible_conda_paths) {
      expanded_path <- path.expand(path)
      if (file.exists(expanded_path)) {
        conda_python <- expanded_path
        break
      }
    }
    
    # If conda environment not found, try to use system python with a warning
    if (is.null(conda_python)) {
      warning("Conda environment 'celltypist_benchmark' not found. Trying system python.")
      conda_python <- "python"
    }
    
    message("Using Python: ", conda_python)
    
    # Construct system command using direct python path
    cmd <- sprintf(
      '"%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"',
      as.character(conda_python),
      as.character(python_script),
      as.character(train_file),
      as.character(test_file), 
      as.character(labels_file),
      as.character(genes_file),
      as.character(output_file),
      as.character(status_file)
    )
    
    message("Executing Python script...")
    message("Command: ", cmd)

    # Execute Python script and capture output
    python_output <- system(cmd, intern = TRUE)
    system_result <- attr(python_output, "status")
    if (is.null(system_result)) system_result <- 0

    # Parse runtime and peak system memory from Python output
    runtime_secs <- NA
    peak_system_memory_mb <- NA
    if (!is.null(python_output)) {
      runtime_line <- grep("RUNTIME_SECS:", python_output, value = TRUE)
      if (length(runtime_line) > 0) {
        runtime_secs <- as.numeric(sub(".*RUNTIME_SECS:([0-9.]+).*", "\\1", runtime_line[1]))
        message(sprintf("Core algorithm runtime: %.4f seconds", runtime_secs))
      }
      memory_line <- grep("PEAK_SYSTEM_MEMORY_MB:", python_output, value = TRUE)
      if (length(memory_line) > 0) {
        peak_system_memory_mb <- as.numeric(sub(".*PEAK_SYSTEM_MEMORY_MB:([0-9.]+).*", "\\1", memory_line[1]))
        message(sprintf("Peak system memory: %.2f MB", peak_system_memory_mb))
      }
    }
    
    # Check execution status
    if (system_result != 0) {
      warning("Python script execution failed with exit code: ", system_result)
      
      # Try to read status file for error details
      if (file.exists(status_file)) {
        status_content <- readLines(status_file, warn = FALSE)
        warning("Python error details: ", paste(status_content, collapse = " "))
      }
      
      return(default_return())
    }
    
    message("✓ Python script completed successfully")
    
    # Check for output files
    if (!file.exists(output_file)) {
      warning("Python output file not found: ", output_file)
      return(default_return())
    }
    
    if (!file.exists(status_file)) {
      warning("Python status file not found: ", status_file)
      return(default_return())
    }
    
    # Read status file
    status_content <- readLines(status_file, warn = FALSE)
    
    if (length(status_content) == 0 || status_content[1] != "SUCCESS") {
      warning("Python script reported failure: ", paste(status_content, collapse = " "))
      return(default_return())
    }
    
    message("✓ Python script reported success")
    
    # Read prediction results
    message("Reading prediction results...")
    
    predictions_df <- read.csv(output_file, stringsAsFactors = FALSE)
    
    # Validate predictions format
    required_columns <- c("cell_id", "predicted_label", "confidence_score")
    missing_columns <- required_columns[!required_columns %in% colnames(predictions_df)]
    
    if (length(missing_columns) > 0) {
      warning("Missing columns in predictions: ", paste(missing_columns, collapse = ", "))
      return(default_return())
    }
    
    # Extract results
    cell_predictions <- predictions_df$predicted_label
    confidence_scores <- predictions_df$confidence_score
    result_cell_ids <- predictions_df$cell_id
    
    # Get true labels from test data
    true_labels <- seurat_test$Ground_Truth_Celltype
    
    # Validate result dimensions
    if (length(cell_predictions) != ncol(seurat_test)) {
      warning(paste("Prediction count mismatch. Expected:", ncol(seurat_test), 
                    "Got:", length(cell_predictions)))
      return(default_return())
    }
    
    # Handle NA values
    cell_predictions[is.na(cell_predictions)] <- "Unknown"
    confidence_scores[is.na(confidence_scores)] <- 0
    
    message("✓ Results processed successfully")
    
    # Clean up temporary files
    temp_files <- c(train_file, test_file, labels_file, genes_file, output_file, status_file)
    
    for (temp_file in temp_files) {
      if (file.exists(temp_file)) {
        tryCatch({
          file.remove(temp_file)
        }, error = function(e) {
          warning("Failed to remove temporary file: ", temp_file)
        })
      }
    }
    
    message("✓ Temporary files cleaned up")
    
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
    warning("CellTypist approach error: ", e$message)
    return(default_return())
  })
}

# For backward compatibility
run_CellTypist <- run_CellTypist_function