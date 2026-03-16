# run_scAnno.R
#################################################
# scAnno Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scAnno Cell Type Annotation Function
#'
#' Purpose: Run scAnno algorithm using deconvolution strategy with TCGA integration
#' Inputs:
#'   - seurat_train: Training Seurat object (used as reference)
#'   - seurat_test: Test Seurat object to predict
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scAnno's deconvolution strategy with TCGA bulk data integration
#'
#' Note: scAnno works at cluster level, so predictions are mapped to individual cells
run_scAnno_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("scAnno", quietly = TRUE)) {
    stop("scAnno package not available. Please install it first.")
  }

  library(scAnno)
  library(Seurat)
  library(dplyr)

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
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
    warning("Ground_Truth_Celltype not found in training data")
    return(default_return())
  }

  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in test data")
    return(default_return())
  }

  tryCatch({


    # Prepare reference expression matrix
    ref.expr <- GetAssayData(seurat_train, assay = "RNA", layer = "data") %>% as.data.frame

    # Set up reference annotations
    Idents(seurat_train) <- "Ground_Truth_Celltype"
    ref.anno <- Idents(seurat_train) %>% as.character

    # Ensure test object has cluster information for scAnno
    # scAnno uses cluster.col to identify clusters in query object
    Idents(seurat_test) <- "Ground_Truth_Celltype"

    # Track peak memory usage
    peak_memory_mb <- NA

    # Run scAnno annotation with memory tracking
    cat("Running scAnno with TCGA integration...\n")
    obj.seu <- seurat_test

    if (!requireNamespace("bench", quietly = TRUE)) {
      warning("bench package not available for memory tracking")
      results <- scAnno(
        query = obj.seu,
        ref.expr = ref.expr,
        ref.anno = ref.anno,
        save.markers = NULL,
        cluster.col = "Ground_Truth_Celltype",
        factor.size = 0.1,
        pvalue.cut = 0.01,
        seed.num = 10,
        redo.markers = FALSE,
        gene.anno = data(gene.anno),
        permut.num = 100,
        permut.p = 0.01,
        show.plot = FALSE,
        verbose = TRUE,
      )
    } else {
      library(bench)
      bench_result <- bench::mark(
        {
          results <- scAnno(
            query = obj.seu,
            ref.expr = ref.expr,
            ref.anno = ref.anno,
            save.markers = NULL,
            cluster.col = "Ground_Truth_Celltype",
            factor.size = 0.1,
            pvalue.cut = 0.01,
            seed.num = 10,
            redo.markers = FALSE,
            gene.anno = data(gene.anno),
            permut.num = 100,
            permut.p = 0.01,
            show.plot = FALSE,
            verbose = TRUE,
          )
        },
        memory = TRUE,
        iterations = 1,
        check = FALSE
      )
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
    }

    # Extract cluster-level results
    truth_clusters <- names(results$pred.label)
    predicted_clusters <- as.vector(results$pred.label)
    cluster_confidence <- results$pred.score

    # Create cluster prediction lookup table
    cluster_prediction_map <- setNames(predicted_clusters, truth_clusters)
    cluster_confidence_map <- setNames(cluster_confidence, truth_clusters)

    # Map cluster predictions to individual cells
    test_cell_clusters <- seurat_test$Ground_Truth_Celltype

    # Get predictions for each cell based on their cluster
    cell_predictions <- cluster_prediction_map[test_cell_clusters]

    # Get confidence scores for each cell based on their cluster
    cell_confidence_scores <- cluster_confidence_map[test_cell_clusters]

    # Handle NA values (clusters not found in prediction map)
    cell_predictions[is.na(cell_predictions)] <- "Unknown"
    cell_confidence_scores[is.na(cell_confidence_scores)] <- 0

    # Get true labels for individual cells
    true_labels <- seurat_test$Ground_Truth_Celltype

    # Validate result dimensions
    if (length(cell_predictions) != ncol(seurat_test)) {
      warning(paste("Prediction count mismatch. Expected:", ncol(seurat_test),
                    "Got:", length(cell_predictions)))
      return(default_return())
    }

    # Return standardized format
    return(list(
      predictions = as.character(cell_predictions),
      true_labels = as.character(true_labels),
      confidence_scores = as.numeric(cell_confidence_scores),
      cell_ids = as.character(colnames(seurat_test)),
      peak_memory_mb = peak_memory_mb
    ))

  }, error = function(e) {
    warning("scAnno error: ", e$message)
    return(default_return())
  })
}

# For backward compatibility
run_scAnno <- run_scAnno_function