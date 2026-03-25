# run_scLearn.R
#################################################
# scLearn Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scLearn Cell Type Annotation Function
#'
#' Purpose: Run scLearn algorithm using discriminative component analysis for learning-based annotation
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference for model training)
#'   - seurat_test: Test Seurat object to predict
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scLearn's DCA-based learning with automatic threshold determination
run_scLearn_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("scLearn", quietly = TRUE)) {
    stop("scLearn package not available. Please install it first: devtools::install_github('bm2-lab/scLearn')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package not available. Please install it first.")
  }
  if (!requireNamespace("M3Drop", quietly = TRUE)) {
    stop("M3Drop package not available. Please install it first.")
  }

  library(scLearn)
  library(SingleCellExperiment)
  library(M3Drop)
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

    cat("Converting Seurat objects to SingleCellExperiment format...\n")

    # Convert Seurat objects to SingleCellExperiment
    sce_train <- as.SingleCellExperiment(seurat_train)
    sce_test <- as.SingleCellExperiment(seurat_test)

    # Extract raw count matrices (scLearn expects raw counts)
    rawcounts_train <- tryCatch({
      # Try to get raw counts first
      GetAssayData(seurat_train, assay = "RNA", layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(seurat_train, assay = "RNA", layer = "counts")
      }, error = function(e2) {
        # Fallback to normalized data if raw counts not available
        GetAssayData(seurat_train, assay = "RNA", layer = "data")
      })
    })

    rawcounts_test <- tryCatch({
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
    if (inherits(rawcounts_train, "sparseMatrix")) {
      rawcounts_train <- as.matrix(rawcounts_train)
    }
    if (inherits(rawcounts_test, "sparseMatrix")) {
      rawcounts_test <- as.matrix(rawcounts_test)
    }

    # Create reference annotations
    refe_ann <- as.character(seurat_train$Ground_Truth_Celltype)
    names(refe_ann) <- colnames(seurat_train)

    # Check for sufficient cells per cell type
    cell_type_counts <- table(refe_ann)
    min_cells_per_type <- 10  # Buffer for Cell_qc() attrition
    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells_per_type]

    if (length(valid_cell_types) < 2) {
      warning("Insufficient cell types with enough cells for scLearn training")
      return(default_return())
    }

    # Filter to keep only valid cell types
    valid_cells <- names(refe_ann)[refe_ann %in% valid_cell_types]
    rawcounts_train_filtered <- rawcounts_train[, valid_cells]
    refe_ann_filtered <- refe_ann[valid_cells]

    cat(sprintf("Training with %d cells across %d cell types\n",
                length(refe_ann_filtered), length(valid_cell_types)))

    # Track peak memory usage
    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Helper to run the core scLearn pipeline
    run_scLearn_pipeline <- function() {
      # Step 1: Cell quality control for training data
      cat("Performing quality control on training data...\n")
      data_qc <- Cell_qc(rawcounts_train_filtered, refe_ann_filtered, species = "Hs",
                         gene_low = 50, umi_low = 50)

      # Step 2: Re-filter post-QC to drop cell types now below threshold
      post_qc_counts <- table(data_qc$sample_information_cellType)
      surviving_types <- names(post_qc_counts)[post_qc_counts >= min_cells_per_type]

      if (length(surviving_types) < 2) {
        warning("Too few cell types remaining after Cell_qc(); returning Unknown")
        return(NULL)
      }

      keep_mask <- data_qc$sample_information_cellType %in% surviving_types
      data_qc$expression_profile          <- data_qc$expression_profile[, keep_mask]
      data_qc$sample_information_cellType <- data_qc$sample_information_cellType[keep_mask]

      # Step 3: Cell type filtering (remove rare cell types)
      cat("Filtering rare cell types...\n")
      data_type_filtered <- Cell_type_filter(
        data_qc$expression_profile,
        data_qc$sample_information_cellType,
        min_cell_number = min_cells_per_type
      )

      # Check if we still have enough cell types after filtering
      remaining_cell_types <- unique(data_type_filtered$sample_information_cellType)
      if (length(remaining_cell_types) < 2) {
        warning("Too few cell types remaining after Cell_type_filter(); returning Unknown")
        return(NULL)
      }

      # Step 4: Feature selection using M3Drop
      cat("Selecting high variance genes...\n")
      high_varGene_names <- tryCatch({
        Feature_selection_M3Drop(data_type_filtered$expression_profile)
      }, error = function(e) {
        warning("M3Drop feature selection failed, using most variable genes")
        gene_vars <- apply(data_type_filtered$expression_profile, 1, var)
        names(gene_vars)[order(gene_vars, decreasing = TRUE)[1:min(2000, length(gene_vars))]]
      })

      cat(sprintf("Selected %d high variance genes\n", length(high_varGene_names)))

      # Step 5: Model training with scLearn
      cat("Training scLearn model...\n")
      scLearn_model_learning_result <- scLearn_model_learning(
        high_varGene_names,
        data_type_filtered$expression_profile,
        data_type_filtered$sample_information_cellType,
        bootstrap_times = 10
      )

      # Step 6: Quality control for test data
      cat("Performing quality control on test data...\n")
      data_qc_query <- Cell_qc(
        rawcounts_test,
        species = "Hs",
        gene_low = 50,
        umi_low = 50
      )

      # Step 7: Cell assignment
      cat("Performing cell type assignment...\n")
      scLearn_cell_assignment(
        scLearn_model_learning_result,
        data_qc_query$expression_profile,
        diff = 0.05,
        threshold_use = TRUE,
        vote_rate = 0.6
      )
    }

    # Run scLearn pipeline with memory tracking
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
      scLearn_predict_result <- tryCatch(run_scLearn_pipeline(), error = function(e) {
        warning(e$message); return(NULL)
      })
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        scLearn_predict_result <- tryCatch(run_scLearn_pipeline(), error = function(e) {
          warning(e$message); return(NULL)
        })
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    if (is.null(scLearn_predict_result)) {
      warning("scLearn pipeline returned NULL — too few cell types after QC/filtering")
      return(default_return())
    }

    # Extract predictions + match
    meta <- seurat_test@meta.data$Ground_Truth_Celltype
    names(meta) <- rownames(seurat_test@meta.data)  # names = barcodes
    
    # Match predicted barcodes to metadata barcodes
    idx <- match(scLearn_predict_result$Query_cell_id, names(meta))
    
    # Get aligned ground truth and predicted labels
    ground_truth <- meta[idx]
    predicted    <- scLearn_predict_result$Predict_cell_type

    # Handle unassigned cells (scLearn might return "unassigned")
    predicted[predicted == "unassigned"] <- "Unknown"
    predicted[is.na(predicted)] <- "Unknown"

    # Create confidence scores
    # scLearn doesn't provide explicit confidence scores, use heuristic
    confidence_scores <- rep(0.8, length(predicted))
    confidence_scores[predicted == "Unknown"] <- 0

    # Return standardized format
    return(list(
      predictions = as.character(predicted),
      true_labels = as.character(ground_truth),
      confidence_scores = as.numeric(confidence_scores),
      cell_ids = as.character(scLearn_predict_result$Query_cell_id[!is.na(idx)]),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    ))

  }, error = function(e) {
    warning("scLearn error: ", e$message)
    return(default_return())
  })
}

# For backward compatibility
run_scLearn <- run_scLearn_function