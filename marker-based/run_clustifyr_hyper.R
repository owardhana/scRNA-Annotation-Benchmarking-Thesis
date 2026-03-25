# run_clustifyr_hyper.R
#################################################
# clustifyr_hyper Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' clustifyr_hyper Cell Type Annotation Function
#'
#' Purpose: Run clustifyr algorithm using hypergeometric enrichment test
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict (will create clusters if needed)
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts markers to ranked list format (like cbmc_m), runs hypergeometric test
#' Note: Uses data frame format with gene names as values (NOT binary matrix)
run_clustifyr_hyper_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("clustifyr", quietly = TRUE)) {
    stop("clustifyr package not available. Please install it with: BiocManager::install('clustifyr')")
  }
  library(clustifyr)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required. Please install it with: install.packages('dplyr')")
  }
  library(dplyr)

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package required. Please install it with: install.packages('tidyr')")
  }
  library(tidyr)

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
    warning("Ground_Truth_Celltype not found in metadata")
    return(default_return())
  }

  if (!is.data.frame(markers) || nrow(markers) == 0) {
    warning("Invalid or empty markers dataframe")
    return(default_return())
  }

  cat("\n=== clustifyr_hyper: Hypergeometric Enrichment Annotation ===\n")

  # Step 1: Ensure test data has clusters
  if (!"seurat_clusters" %in% colnames(seurat_test@meta.data)) {
    cat("Creating clusters for test data...\n")
    tryCatch({
      seurat_test <- FindNeighbors(seurat_test, verbose = FALSE)
      seurat_test <- FindClusters(seurat_test, resolution = 0.5, verbose = FALSE)
    }, error = function(e) {
      warning(paste("Failed to create clusters:", e$message))
      return(default_return())
    })
  }

  n_test_clusters <- length(unique(seurat_test$seurat_clusters))
  cat(sprintf("Test data: %d clusters across %d cells\n",
              n_test_clusters, ncol(seurat_test)))

  # Step 2: Convert markers to ranked list format (like cbmc_m)
  # Format: data frame with rows = ranks, columns = cell types, values = gene names
  cat("\nConverting markers to ranked list format (data frame with gene names)...\n")

  # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
  filtered_markers <- markers[markers$avg_log2FC >= 0.5 &
                              markers$p_val_adj < 0.05 &
                              markers$pct.1 >= 0.15, ]

  if (nrow(filtered_markers) == 0) {
    warning("No high-quality markers found")
    return(default_return())
  }

  cat(sprintf("Filtered markers: %d genes across %d cell types\n",
              nrow(filtered_markers),
              length(unique(filtered_markers$cluster))))

  # Take top 50 markers per cell type (standardized)
  markers_per_type <- 50
  top_markers <- filtered_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = markers_per_type) %>%
    dplyr::mutate(rank = row_number()) %>%
    dplyr::ungroup()

  cat(sprintf("Using top %d markers per cell type (standardized)\n", markers_per_type))
  cat("Markers per cell type:\n")
  print(table(top_markers$cluster))

  # Convert to wide format: rows = ranks, columns = cell types, values = gene names
  # This matches the cbmc_m format
  marker_list <- top_markers %>%
    dplyr::select(cluster, gene, rank) %>%
    tidyr::pivot_wider(
      names_from = cluster,
      values_from = gene,
      values_fill = NA
    ) %>%
    dplyr::select(-rank) %>%
    as.data.frame()

  # Remove row names (just use numeric indices)
  rownames(marker_list) <- NULL

  cat(sprintf("\nRanked marker list created: %d ranks x %d cell types\n",
              nrow(marker_list), ncol(marker_list)))
  cat("Format check:\n")
  cat(sprintf("  Class: %s\n", class(marker_list)))
  cat(sprintf("  Column 1 type: %s\n", class(marker_list[,1])))
  cat("\nSample marker list (first 10 rows):\n")
  print(head(marker_list, 10))

  # Step 3: Run clustifyr with hypergeometric test and memory tracking
  cat("\nRunning clustify_lists with hypergeometric test...\n")

  runtime_secs <- NA
  peak_system_memory_mb <- NA
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()
    seurat_annotated <- tryCatch({
      clustify_lists(
        input = seurat_test,
        marker = marker_list,
        cluster_col = "seurat_clusters",
        metric = "hyper",
        obj_out = TRUE
      )
    }, error = function(e) {
      warning(paste("clustify_lists failed:", e$message))
      cat("Error details:", e$message, "\n")
      return(NULL)
    })
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        seurat_annotated <- clustify_lists(
          input = seurat_test,
          marker = marker_list,
          cluster_col = "seurat_clusters",
          metric = "hyper",
          obj_out = TRUE
        )
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  if (is.null(seurat_annotated)) {
    return(default_return())
  }

  # Step 4: Extract predictions from 'type' column
  if (!"type" %in% colnames(seurat_annotated@meta.data)) {
    warning("clustifyr did not create 'type' column in metadata")
    return(default_return())
  }

  predictions <- as.character(seurat_annotated$type)
  true_labels <- seurat_annotated$Ground_Truth_Celltype

  # Step 5: Get confidence scores
  # clustifyr doesn't return confidence scores by default
  # We'll set all to 1.0 for successful predictions, 0.0 for Unknown
  confidence_scores <- ifelse(predictions == "Unknown", 0, 1)

  # Summary
  cat(sprintf("\nPrediction summary:\n"))
  cat(sprintf("  Total cells: %d\n", length(predictions)))
  cat(sprintf("  Assigned: %d\n", sum(predictions != "Unknown")))
  cat(sprintf("  Unknown: %d\n", sum(predictions == "Unknown")))
  cat(sprintf("  Unique predicted types: %d\n",
              length(unique(predictions[predictions != "Unknown"]))))

  cat("\nPredictions distribution:\n")
  print(table(predictions))

  # Return standardized format
  return(list(
    predictions = as.character(predictions),
    true_labels = true_labels,
    confidence_scores = as.numeric(confidence_scores),
    cell_ids = colnames(seurat_annotated),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
}

# For backward compatibility
run_clustifyr_hyper <- run_clustifyr_hyper_function
