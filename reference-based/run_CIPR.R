# run_CIPR.R
#################################################
# CIPR Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' CIPR (Cluster Identity Predictor) Cell Type Annotation Function
#'
#' Purpose: Run CIPR algorithm using logFC dot product comparison between clusters
#' Inputs:
#'   - seurat_train: Training Seurat object (used as custom reference)
#'   - seurat_test: Test Seurat object to predict (will be clustered if needed)
#'   - markers: Marker genes dataframe from FindAllMarkers() (not used - recomputed internally)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: CIPR is cluster-based. Steps:
#'   1. Cluster test data (if not already clustered)
#'   2. Find markers for test clusters
#'   3. Build custom reference from training data average expression
#'   4. Run CIPR with logfc_dot_product method
#'   5. Map cluster predictions to individual cells
#' Note: CIPR creates global objects that are cleaned up after extraction
run_CIPR_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("CIPR", quietly = TRUE)) {
    stop("CIPR package not available. Please install it with: devtools::install_github('atakanekiz/CIPR-Package')")
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required. Install with: install.packages('dplyr')")
  }

  library(CIPR)
  library(dplyr)

  # Default return for errors
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

  # Validate inputs
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in test metadata")
    return(default_return())
  }

  if (!"Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
    warning("Ground_Truth_Celltype not found in train metadata")
    return(default_return())
  }

  cat("=== CIPR Reference-Based Annotation ===\n")

  # Step 1: Use Ground_Truth_Celltype as clusters
  # CIPR is cluster-based - we use the true cell types as cluster assignments
  # This allows CIPR to find markers for each cell type group and compare to reference
  cat("Using Ground_Truth_Celltype as cluster assignments...\n")

  n_test_clusters <- length(unique(seurat_test$Ground_Truth_Celltype))
  cat(sprintf("Test data: %d cell type groups across %d cells\n",
              n_test_clusters, ncol(seurat_test)))

  # Track peak memory usage
  runtime_secs <- NA
  peak_system_memory_mb <- NA

  # Run CIPR pipeline with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
    start_time <- Sys.time()

    # Step 2: Find markers for test cell type groups
    cat("Finding markers for test cell type groups...\n")
    original_idents <- Idents(seurat_test)
    Idents(seurat_test) <- seurat_test$Ground_Truth_Celltype

    test_markers <- tryCatch({
      FindAllMarkers(seurat_test, only.pos = TRUE, min.pct = 0.1,
                     logfc.threshold = 0.1, verbose = FALSE)
    }, error = function(e) {
      warning(paste("FindAllMarkers failed for test data:", e$message))
      return(NULL)
    })

    Idents(seurat_test) <- original_idents

    if (is.null(test_markers) || nrow(test_markers) == 0) {
      warning("No markers found for test clusters")
      return(default_return())
    }

    cat(sprintf("Found %d marker genes across %d clusters\n",
                nrow(test_markers), length(unique(test_markers$cluster))))

    # Step 3: Build custom reference
    cat("Building custom reference from training data...\n")
    custom_ref_expr <- tryCatch({
      avgexp <- AverageExpression(seurat_train, assays = "RNA", verbose = FALSE)$RNA
      ref_df <- as.data.frame(avgexp)
      ref_df$gene <- tolower(rownames(avgexp))
      ref_df <- ref_df[, c("gene", setdiff(colnames(ref_df), "gene"))]
      ref_df
    }, error = function(e) {
      warning(paste("Failed to create reference expression matrix:", e$message))
      return(NULL)
    })

    if (is.null(custom_ref_expr) || nrow(custom_ref_expr) == 0) {
      warning("Reference expression matrix is empty")
      return(default_return())
    }

    # Step 4: Build reference metadata
    train_meta <- seurat_train@meta.data
    if (!"seurat_clusters" %in% colnames(train_meta)) {
      train_meta$seurat_clusters <- train_meta$Ground_Truth_Celltype
    }

    cluster_celltype_map <- train_meta %>%
      dplyr::select(seurat_clusters, Ground_Truth_Celltype) %>%
      distinct() %>%
      arrange(seurat_clusters)

    ref_columns <- setdiff(colnames(custom_ref_expr), "gene")
    custom_ref_annot <- data.frame(short_name = ref_columns, stringsAsFactors = FALSE)

    if (all(ref_columns %in% cluster_celltype_map$Ground_Truth_Celltype)) {
      custom_ref_annot$long_name <- custom_ref_annot$short_name
      custom_ref_annot$reference_cell_type <- custom_ref_annot$short_name
    } else {
      custom_ref_annot$long_name <- cluster_celltype_map$Ground_Truth_Celltype[
        match(custom_ref_annot$short_name, as.character(cluster_celltype_map$seurat_clusters))
      ]
      custom_ref_annot$reference_cell_type <- custom_ref_annot$long_name
    }

    custom_ref_annot$description <- "Training reference"
    custom_ref_annot$long_name[is.na(custom_ref_annot$long_name)] <- custom_ref_annot$short_name[is.na(custom_ref_annot$long_name)]
    custom_ref_annot$reference_cell_type[is.na(custom_ref_annot$reference_cell_type)] <- custom_ref_annot$short_name[is.na(custom_ref_annot$reference_cell_type)]

    cat(sprintf("Reference: %d cell types with %d genes\n",
                nrow(custom_ref_annot), nrow(custom_ref_expr)))

    # Step 5: Run CIPR
    cat("Running CIPR with logfc_dot_product method...\n")
    cipr_success <- tryCatch({
      CIPR(input_dat = test_markers, comp_method = "logfc_dot_product", reference = "custom",
           custom_reference = custom_ref_expr, custom_ref_annot = custom_ref_annot,
           keep_top_var = 100, plot_ind = FALSE, plot_top = FALSE,
           global_results_obj = TRUE, global_plot_obj = FALSE)
      TRUE
    }, error = function(e) {
      warning(paste("CIPR execution failed:", e$message))
      FALSE
    })

    if (!cipr_success) {
      return(default_return())
    }

    # Step 6: Extract results
    if (!exists("CIPR_all_results", envir = .GlobalEnv)) {
      warning("CIPR did not create results object")
      return(default_return())
    }

    cipr_results <- get("CIPR_all_results", envir = .GlobalEnv)
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)

    run_cipr_pipeline <- function() {
      # Step 2: Find markers for test cell type groups
      cat("Finding markers for test cell type groups...\n")
      original_idents <- Idents(seurat_test)
      Idents(seurat_test) <- seurat_test$Ground_Truth_Celltype

      test_markers <- tryCatch({
        FindAllMarkers(seurat_test, only.pos = TRUE, min.pct = 0.1,
                       logfc.threshold = 0.1, verbose = FALSE)
      }, error = function(e) {
        warning(paste("FindAllMarkers failed for test data:", e$message))
        return(NULL)
      })

      Idents(seurat_test) <- original_idents

      if (is.null(test_markers) || nrow(test_markers) == 0) {
        stop("No markers found for test clusters")
      }

      cat(sprintf("Found %d marker genes across %d clusters\n",
                  nrow(test_markers), length(unique(test_markers$cluster))))

      # Step 3: Build custom reference
      cat("Building custom reference from training data...\n")
      custom_ref_expr <- tryCatch({
        avgexp <- AverageExpression(seurat_train, assays = "RNA", verbose = FALSE)$RNA
        ref_df <- as.data.frame(avgexp)
        ref_df$gene <- tolower(rownames(avgexp))
        ref_df <- ref_df[, c("gene", setdiff(colnames(ref_df), "gene"))]
        ref_df
      }, error = function(e) {
        warning(paste("Failed to create reference expression matrix:", e$message))
        return(NULL)
      })

      if (is.null(custom_ref_expr) || nrow(custom_ref_expr) == 0) {
        stop("Reference expression matrix is empty")
      }

      # Step 4: Build reference metadata
      train_meta <- seurat_train@meta.data
      if (!"seurat_clusters" %in% colnames(train_meta)) {
        train_meta$seurat_clusters <- train_meta$Ground_Truth_Celltype
      }

      cluster_celltype_map <- train_meta %>%
        dplyr::select(seurat_clusters, Ground_Truth_Celltype) %>%
        distinct() %>%
        arrange(seurat_clusters)

      ref_columns <- setdiff(colnames(custom_ref_expr), "gene")
      custom_ref_annot <- data.frame(short_name = ref_columns, stringsAsFactors = FALSE)

      if (all(ref_columns %in% cluster_celltype_map$Ground_Truth_Celltype)) {
        custom_ref_annot$long_name <- custom_ref_annot$short_name
        custom_ref_annot$reference_cell_type <- custom_ref_annot$short_name
      } else {
        custom_ref_annot$long_name <- cluster_celltype_map$Ground_Truth_Celltype[
          match(custom_ref_annot$short_name, as.character(cluster_celltype_map$seurat_clusters))
        ]
        custom_ref_annot$reference_cell_type <- custom_ref_annot$long_name
      }

      custom_ref_annot$description <- "Training reference"
      custom_ref_annot$long_name[is.na(custom_ref_annot$long_name)] <- custom_ref_annot$short_name[is.na(custom_ref_annot$long_name)]
      custom_ref_annot$reference_cell_type[is.na(custom_ref_annot$reference_cell_type)] <- custom_ref_annot$short_name[is.na(custom_ref_annot$reference_cell_type)]

      cat(sprintf("Reference: %d cell types with %d genes\n",
                  nrow(custom_ref_annot), nrow(custom_ref_expr)))

      # Step 5: Run CIPR
      cat("Running CIPR with logfc_dot_product method...\n")
      cipr_success <- tryCatch({
        CIPR(input_dat = test_markers, comp_method = "logfc_dot_product", reference = "custom",
             custom_reference = custom_ref_expr, custom_ref_annot = custom_ref_annot,
             keep_top_var = 100, plot_ind = FALSE, plot_top = FALSE,
             global_results_obj = TRUE, global_plot_obj = FALSE)
        TRUE
      }, error = function(e) {
        warning(paste("CIPR execution failed:", e$message))
        FALSE
      })

      if (!cipr_success) {
        stop("CIPR execution failed")
      }

      # Step 6: Extract results
      if (!exists("CIPR_all_results", envir = .GlobalEnv)) {
        stop("CIPR did not create results object")
      }

      get("CIPR_all_results", envir = .GlobalEnv)
    }

    peakRAM_result <- peakRAM::peakRAM({
      cipr_results <- tryCatch(run_cipr_pipeline(), error = function(e) {
        warning(e$message); return(NULL)
      })
    })

    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]

    if (is.null(cipr_results)) {
      warning("CIPR pipeline failed during memory tracking")
      return(default_return())
    }
  }

  if (is.null(cipr_results) || nrow(cipr_results) == 0) {
    warning("CIPR results are empty")
    return(default_return())
  }

  cat(sprintf("CIPR completed: %d cluster-reference comparisons\n", nrow(cipr_results)))

  # Step 7: Get top prediction per cluster
  cluster_predictions <- tryCatch({
    cipr_results %>%
      group_by(cluster) %>%
      slice_max(order_by = identity_score, n = 1, with_ties = FALSE) %>%
      dplyr::select(cluster, predicted_type = reference_cell_type,
             identity_score, z_score) %>%
      ungroup()
  }, error = function(e) {
    warning(paste("Failed to extract cluster predictions:", e$message))
    return(NULL)
  })

  if (is.null(cluster_predictions) || nrow(cluster_predictions) == 0) {
    warning("Could not extract cluster predictions from CIPR results")
    return(default_return())
  }

  # Normalize predicted_type to match original Ground_Truth_Celltype names
  # CIPR may internally transform names (e.g., underscores to hyphens)
  original_types <- unique(as.character(seurat_train$Ground_Truth_Celltype))
  normalize_name <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
  name_lookup <- setNames(original_types, normalize_name(original_types))
  normalized_preds <- normalize_name(as.character(cluster_predictions$predicted_type))
  cluster_predictions$predicted_type <- ifelse(
    normalized_preds %in% names(name_lookup),
    name_lookup[normalized_preds],
    as.character(cluster_predictions$predicted_type)
  )

  cat("Top predictions per cell type group:\n")
  print(as.data.frame(cluster_predictions))

  # Step 8: Map cluster predictions to individual cells
  # Use Ground_Truth_Celltype as the grouping (since we used it for markers)
  test_cell_types <- as.character(seurat_test$Ground_Truth_Celltype)

  # Match each cell's true cell type to the cluster prediction
  cell_predictions <- cluster_predictions$predicted_type[
    match(test_cell_types, as.character(cluster_predictions$cluster))
  ]

  # Get confidence scores (z-scores)
  confidence_scores <- cluster_predictions$z_score[
    match(test_cell_types, as.character(cluster_predictions$cluster))
  ]

  # Step 9: Handle unmatched clusters
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0

  # Step 10: Clean up global environment
  # CIPR pollutes the global environment - clean up
  if (exists("CIPR_all_results", envir = .GlobalEnv)) {
    rm(CIPR_all_results, envir = .GlobalEnv)
  }
  if (exists("CIPR_top_results", envir = .GlobalEnv)) {
    rm(CIPR_top_results, envir = .GlobalEnv)
  }

  # Get true labels
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Summary statistics
  cat(sprintf("\nPrediction summary:\n"))
  cat(sprintf("  Total cells: %d\n", length(cell_predictions)))
  cat(sprintf("  Assigned: %d\n", sum(cell_predictions != "Unknown")))
  cat(sprintf("  Unknown: %d\n", sum(cell_predictions == "Unknown")))
  cat(sprintf("  Unique predicted types: %d\n", length(unique(cell_predictions[cell_predictions != "Unknown"]))))

  cat("\nNote: CIPR used Ground_Truth_Celltype as cluster grouping\n")
  cat("Each cell type group's markers were compared to training reference\n")

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = as.numeric(confidence_scores),
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
}

# For backward compatibility
run_CIPR <- run_CIPR_function
