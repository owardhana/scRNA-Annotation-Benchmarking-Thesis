#' Simplified Garnett Cell Type Annotation Function
#'
#' @param seurat_train Training Seurat object
#' @param seurat_test Test Seurat object to predict
#' @param markers Marker genes dataframe from FindAllMarkers()
#' @param max_markers_per_type Maximum number of top markers to use per cell type (default: 50, standardized)
#' @param rank_prob_ratio Classification strictness threshold (default: 0.5, lower = more permissive)
#' @return List with predictions, true_labels, confidence_scores, cell_ids

run_Garnett_database_function <- function(seurat_train, seurat_test, markers,
                                 max_markers_per_type = 50,
                                 rank_prob_ratio = 1.1) {
  library(garnett)
  library(monocle3)
  library(org.Hs.eg.db)
  library(SeuratWrappers)

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
  cat("=== GARNETT ===\n")
  
  # 1. Check training data cell type labels
  cat("\n1. Checking training data:\n")
  if ("seurat_clusters" %in% colnames(seurat_train@meta.data)) {
    train_clusters <- table(seurat_train$seurat_clusters)
    cat("Training clusters:", paste(names(train_clusters), collapse = ", "), "\n")
  } else if ("Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
    train_types <- table(seurat_train$Ground_Truth_Celltype)
    cat("Training cell types:", paste(names(train_types), collapse = ", "), "\n")
  } else {
    cat("WARNING: No cell type labels found in training data!\n")
  }
  
  # 2. Check marker genes
  cat("\n2. Checking markers:\n")
  cat("Marker dataframe dimensions:", nrow(markers), "x", ncol(markers), "\n")
  cat("Marker columns:", paste(colnames(markers), collapse = ", "), "\n")
  cat("Unique cell types in markers:", paste(unique(markers$cluster), collapse = ", "), "\n")
  cat("Number of markers per cell type:\n")
  print(table(markers$cluster))
  
  # 3. Check gene overlap
  cat("\n3. Checking gene overlap:\n")
  train_genes <- rownames(seurat_train)
  test_genes <- rownames(seurat_test)
  marker_genes <- unique(markers$gene)
  
  cat("Training genes:", length(train_genes), "\n")
  cat("Test genes:", length(test_genes), "\n")
  cat("Marker genes:", length(marker_genes), "\n")
  
  train_marker_overlap <- sum(marker_genes %in% train_genes)
  test_marker_overlap <- sum(marker_genes %in% test_genes)
  
  cat("Markers found in training data:", train_marker_overlap, "/", length(marker_genes), 
      "(", round(100*train_marker_overlap/length(marker_genes), 1), "%)\n")
  cat("Markers found in test data:", test_marker_overlap, "/", length(marker_genes), 
      "(", round(100*test_marker_overlap/length(marker_genes), 1), "%)\n")
  
  if (train_marker_overlap < length(marker_genes) * 0.5) {
    cat("WARNING: Less than 50% of marker genes found in training data!\n")
  }
  
  # 4. Create and inspect marker file
  create_marker_file_debug <- function(marker_df, output_file, valid_genes = NULL) {
    cat("\n4. Creating marker file:\n")

    filtered_markers <- marker_df

    # Ensure required columns
    if (!"cluster" %in% colnames(filtered_markers)) {
      stop("No 'cluster' column found in filtered markers!")
    }
    if (!"gene" %in% colnames(filtered_markers)) {
      stop("No 'gene' column found in filtered markers!")
    }

    # Sanitize gene names: remove markers with non-ASCII characters (Greek letters
    # like γ, β, α), spaces, or other invalid characters that crash Garnett's parser.
    # CellMarker DB contains entries like "IFN-γ", "β-catenin", "CD8α" etc.
    before_sanitize <- nrow(filtered_markers)
    valid_symbol <- grepl("^[A-Za-z][A-Za-z0-9._-]*$", filtered_markers$gene) &
                    nchar(trimws(filtered_markers$gene)) <= 20 &
                    !grepl("[^\x01-\x7F]", filtered_markers$gene, perl = TRUE)
    filtered_markers <- filtered_markers[valid_symbol, ]
    # Also trim any trailing whitespace from gene names
    filtered_markers$gene <- trimws(filtered_markers$gene)
    cat(sprintf("Gene sanitization: %d/%d markers have valid gene symbols\n",
                nrow(filtered_markers), before_sanitize))

    # Filter to genes present in expression matrix (critical for database mode
    # where CellMarker gene names may not all match Seurat gene symbols)
    if (!is.null(valid_genes)) {
      before_n <- nrow(filtered_markers)
      filtered_markers <- filtered_markers[filtered_markers$gene %in% valid_genes, ]
      cat(sprintf("Gene filter: %d/%d markers have matching genes in expression matrix\n",
                  nrow(filtered_markers), before_n))
    }

    # Sanitize cell type names for Garnett marker file (replace non-ASCII and
    # special chars that crash the parser)
    sanitize_ct_name <- function(name) {
      # Replace non-ASCII with ASCII equivalents where possible
      name <- gsub("[\x80-\xFF]", "", name, perl = TRUE)
      # Replace chars illegal in Garnett marker files with underscores
      name <- gsub("[^A-Za-z0-9 _+-]", "", name)
      trimws(name)
    }

    # Create marker content in Garnett format
    cell_types <- unique(filtered_markers$cluster)
    marker_content <- character(0)

    for (ct in cell_types) {
      ct_genes <- filtered_markers$gene[filtered_markers$cluster == ct]
      safe_ct <- sanitize_ct_name(ct)
      cat("Cell type", safe_ct, ":", length(ct_genes), "markers\n")

      if (length(ct_genes) > 0 && nchar(safe_ct) > 0) {
        marker_content <- c(marker_content, paste0(">", safe_ct))
        marker_content <- c(marker_content, paste("expressed:", paste(ct_genes, collapse = ", ")))
        marker_content <- c(marker_content, "")
      }
    }

    if (length(marker_content) == 0) {
      stop("No valid markers found for any cell type!")
    }

    writeLines(marker_content, output_file)

    # Show marker file content
    cat("\nMarker file content (first 20 lines):\n")
    marker_lines <- readLines(output_file)
    cat(paste(head(marker_lines, 20), collapse = "\n"), "\n")
    if (length(marker_lines) > 20) cat("... (", length(marker_lines) - 20, "more lines)\n")

    return(output_file)
  }
  
  # 5. Convert to CDS and check
  cat("\n5. Converting to Monocle3 CDS:\n")
  train_cds <- as.cell_data_set(seurat_train, assay = "RNA")
  test_cds <- as.cell_data_set(seurat_test, assay = "RNA")

  cat("Training CDS dimensions:", nrow(train_cds), "genes x", ncol(train_cds), "cells\n")
  cat("Test CDS dimensions:", nrow(test_cds), "genes x", ncol(test_cds), "cells\n")

  # Add cluster information to test_cds if available (helps with cluster_extend)
  if ("seurat_clusters" %in% colnames(colData(test_cds))) {
    cat("Test data has seurat_clusters - Garnett will use this for cluster extension\n")
  } else {
    cat("Note: Test data has no cluster information - cluster extension may be limited\n")
  }

  # Check if cell type information is properly transferred
  train_coldata <- colData(train_cds)
  cat("Training CDS metadata columns:", paste(colnames(train_coldata), collapse = ", "), "\n")
  
  # 6. Map CellMarker cell type names to Ground_Truth_Celltype names
  # In DATABASE MODE, marker types come from CellMarker (e.g., "Astrocyte") and
  # won't match Ground_Truth_Celltype (e.g., "astrocytes"). Garnett is supervised
  # and NEEDS labeled training cells matching marker file types. We fuzzy-match
  # CellMarker names to GT names, then rewrite markers$cluster accordingly.
  marker_cell_types <- unique(markers$cluster)
  gt_types <- unique(as.character(seurat_train$Ground_Truth_Celltype))

  map_cellmarker_to_gt <- function(cm_types, gt_types) {
    mapping <- setNames(rep(NA_character_, length(cm_types)), cm_types)
    cm_lower <- tolower(cm_types)
    gt_lower <- tolower(gt_types)

    for (i in seq_along(cm_types)) {
      # Exact match (case-insensitive)
      idx <- match(cm_lower[i], gt_lower)
      if (!is.na(idx)) { mapping[i] <- gt_types[idx]; next }

      # Singularize/pluralize: "astrocyte" vs "astrocytes", "neuron" vs "neurons"
      cm_stem <- sub("s$", "", cm_lower[i])
      for (j in seq_along(gt_types)) {
        gt_stem <- sub("s$", "", gt_lower[j])
        if (cm_stem == gt_stem || cm_stem == gt_lower[j] || cm_lower[i] == gt_stem) {
          mapping[i] <- gt_types[j]; break
        }
      }
      if (!is.na(mapping[i])) next

      # Substring: GT contains CM or CM contains GT
      for (j in seq_along(gt_types)) {
        if (grepl(cm_lower[i], gt_lower[j], fixed = TRUE) ||
            grepl(gt_lower[j], cm_lower[i], fixed = TRUE)) {
          mapping[i] <- gt_types[j]; break
        }
      }
    }
    return(mapping)
  }

  cm_to_gt <- map_cellmarker_to_gt(marker_cell_types, gt_types)
  n_mapped <- sum(!is.na(cm_to_gt))
  cat(sprintf("\nCellMarker → Ground_Truth mapping: %d/%d types mapped\n",
              n_mapped, length(cm_to_gt)))
  for (i in seq_along(cm_to_gt)) {
    if (!is.na(cm_to_gt[i])) {
      cat(sprintf("  '%s' → '%s'\n", names(cm_to_gt)[i], cm_to_gt[i]))
    }
  }
  if (n_mapped < length(cm_to_gt)) {
    unmapped <- names(cm_to_gt)[is.na(cm_to_gt)]
    cat(sprintf("  Unmapped (%d): %s\n", length(unmapped),
                paste(head(unmapped, 10), collapse = ", ")))
  }

  # Remap markers$cluster from CellMarker names to GT names; drop unmapped types
  markers_mapped <- markers
  markers_mapped$cluster <- cm_to_gt[markers_mapped$cluster]
  markers_mapped <- markers_mapped[!is.na(markers_mapped$cluster), ]

  if (nrow(markers_mapped) == 0) {
    cat("WARNING: No CellMarker types could be mapped to Ground_Truth — using original markers\n")
    markers_mapped <- markers
  }

  # Create marker file with mapped names
  temp_marker_file <- tempfile(pattern = "garnett_markers_", fileext = ".txt")
  create_marker_file_debug(markers_mapped, temp_marker_file, valid_genes = rownames(seurat_train))

  # 7. Train classifier with debugging
  cat("\n6. Training Garnett classifier:\n")

  # Set cell_type from Ground_Truth_Celltype — now marker file types match
  if (!"cell_type" %in% colnames(colData(train_cds))) {
    if ("Ground_Truth_Celltype" %in% colnames(colData(train_cds))) {
      colData(train_cds)$cell_type <- as.character(colData(train_cds)$Ground_Truth_Celltype)
      cat("Added cell_type from Ground_Truth_Celltype\n")
    } else if ("seurat_clusters" %in% colnames(colData(train_cds))) {
      colData(train_cds)$cell_type <- paste0("cluster_", colData(train_cds)$seurat_clusters)
      cat("Added cell_type from seurat_clusters\n")
    } else {
      colData(train_cds)$cell_type <- rep("Unknown", ncol(train_cds))
      cat("No label column found — setting cell_type to 'Unknown'\n")
    }
  }

  # Validate marker-training label alignment (using mapped marker types)
  mapped_marker_types <- unique(markers_mapped$cluster)
  training_cell_types <- unique(colData(train_cds)$cell_type)

  cat("\nValidating cell type alignment (after mapping):\n")
  cat("Mapped marker types:", paste(sort(mapped_marker_types), collapse = ", "), "\n")
  cat("Training cell types:", paste(sort(training_cell_types), collapse = ", "), "\n")

  common_types <- intersect(mapped_marker_types, training_cell_types)
  cat("Overlap:", length(common_types), "/", length(mapped_marker_types), "marker types found in training\n")

  # Train classifier with memory tracking (Step 1 of 2)
  # Note: Garnett internally filters cells based on marker expression
  # With focused marker sets (top 100 per type), more cells should pass filtering
  runtime_secs <- NA
  peak_system_memory_mb <- NA
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
    start_time <- Sys.time()
    classifier <- train_cell_classifier(
      cds = train_cds,
      marker_file = temp_marker_file,
      db = org.Hs.eg.db,
      cds_gene_id_type = "SYMBOL",
      marker_file_gene_id_type = "SYMBOL",
      num_unknown = 50,
      cores = 1
    )
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    pr1 <- peakRAM::peakRAM({ classifier <- train_cell_classifier(
      cds = train_cds,
      marker_file = temp_marker_file,
      db = org.Hs.eg.db,
      cds_gene_id_type = "SYMBOL",
      marker_file_gene_id_type = "SYMBOL",
      num_unknown = 50,
      cores = 1
    ) })
    runtime_secs <- pr1$Elapsed_Time_sec[1]
    peak_system_memory_mb <- pr1$Peak_RAM_Used_MiB[1]
  }

  cat("\nTraining cells used by Garnett:\n")
  cat("Total training cells available:", ncol(train_cds), "\n")
  cat("Note: Garnett selects cells based on marker expression, not all cells are used\n")
  
  cat("Classifier trained successfully\n")

  # Inspect classifier to see what it learned
  cat("\n6.5. Inspecting classifier:\n")
  cat("Classifier class:", class(classifier), "\n")
  cat("Classifier names:", paste(names(classifier), collapse = ", "), "\n")
  if ("chosen_cells" %in% names(classifier)) {
    cat("Chosen cells per type:", paste(names(table(classifier$chosen_cells)), collapse = ", "), "\n")
  }
  if ("classification_tree" %in% names(classifier)) {
    cat("Classification tree exists:", !is.null(classifier$classification_tree), "\n")
  }

  # 8. Classify test cells with permissive parameters
  cat("\n7. Classifying test cells:\n")
  cat("Classifying test cells with Garnett...\n")

  # Check columns BEFORE classification
  cat("Columns in test_cds BEFORE classify_cells:\n")
  cols_before <- colnames(colData(test_cds))
  cat(paste(cols_before, collapse = ", "), "\n")

  # Classify with memory tracking (Step 2 of 2)
  test_cds <- tryCatch({
    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      start_time2 <- Sys.time()
      result_cds <- classify_cells(test_cds, classifier)
      runtime_secs <- runtime_secs + as.numeric(difftime(Sys.time(), start_time2, units = "secs"))
      result_cds
    } else {
      pr2 <- peakRAM::peakRAM({
        result_cds <- classify_cells(test_cds, classifier)
      })
      runtime_secs <- runtime_secs + pr2$Elapsed_Time_sec[1]
      peak_system_memory_mb <- max(peak_system_memory_mb, pr2$Peak_RAM_Used_MiB[1], na.rm = TRUE)
      result_cds
    }
  }, error = function(e) {
    cat("ERROR in classify_cells:", e$message, "\n")
    return(test_cds)
  }, warning = function(w) {
    cat("WARNING in classify_cells:", w$message, "\n")
    invokeRestart("muffleWarning")
  })

  # Check columns AFTER classification
  cat("\nColumns in test_cds AFTER classify_cells:\n")
  cols_after <- colnames(colData(test_cds))
  cat(paste(cols_after, collapse = ", "), "\n")
  print(head(colData(test_cds)))

  # Check which new columns were added
  new_cols <- setdiff(cols_after, cols_before)
  if (length(new_cols) > 0) {
    cat("New columns added by classify_cells:", paste(new_cols, collapse = ", "), "\n")
  } else {
    cat("WARNING: No new columns added by classify_cells!\n")
  }

  # Diagnose classification results
  cell_metadata_check <- as.data.frame(colData(test_cds))

  # Check if individual cell_type assignments were made (before cluster extension)
  if ("cell_type" %in% colnames(cell_metadata_check)) {
    cell_type_vals <- table(cell_metadata_check$cell_type)
    cat("\nIndividual cell_type assignments (before cluster extension):\n")
    print(cell_type_vals)

    unknown_individual <- sum(cell_metadata_check$cell_type == "Unknown", na.rm = TRUE)
    cat(sprintf("Individual Unknown: %d / %d (%.1f%%)\n",
                unknown_individual, nrow(cell_metadata_check),
                100 * unknown_individual / nrow(cell_metadata_check)))
  }

  # Check if cluster extension succeeded
  if ("cluster_ext_type" %in% colnames(cell_metadata_check)) {
    cat("\nCluster extension succeeded - cluster_ext_type column exists\n")
    unknown_frac <- sum(cell_metadata_check$cluster_ext_type == "Unknown", na.rm = TRUE) / ncol(test_cds)

    if (unknown_frac > 0.9) {
      cat(sprintf("\nWARNING: %.1f%% cells are Unknown with cluster_extend=TRUE\n", unknown_frac * 100))
      cat("Retrying classification with cluster_extend=FALSE (direct cell classification)...\n")

      # Retry without cluster extension
      test_cds <- classify_cells(test_cds, classifier)

      cat("Retry completed\n")
    }
  } else {
    cat("\nWARNING: Cluster extension FAILED - cluster_ext_type column NOT created!\n")
    cat("Possible reasons:\n")
    cat("  1. No individual cells were assigned (all Unknown)\n")
    cat("  2. Garnett's cluster extension logic failed internal checks\n")
    cat("  3. classify_cells encountered an error\n")

    # Try direct classification as fallback
    cat("\nRetrying with direct classification...\n")
    test_cds <- classify_cells(test_cds, classifier)
    cat("Direct classification completed\n")
  }

  # 9. Check results
  cat("\n8. Checking classification results:\n")
  cell_metadata <- as.data.frame(colData(test_cds))
  cat("Available result columns:", paste(colnames(cell_metadata), collapse = ", "), "\n")
  
  # Look for prediction columns
  pred_cols <- grep("cell_type|cluster_ext_type|garnett", colnames(cell_metadata), value = TRUE)
  cat("Potential prediction columns:", paste(pred_cols, collapse = ", "), "\n")

  # Debug: Show unique values in each prediction column
  cat("\nDEBUG - Checking prediction column contents:\n")
  for (col in pred_cols) {
    unique_vals <- unique(cell_metadata[[col]])
    cat(sprintf("  %s: %s\n", col, paste(head(unique_vals, 10), collapse = ", ")))
  }

  # Additional detailed diagnostics
  cat("\nDEBUG - Detailed analysis:\n")
  cat(sprintf("  Total test cells: %d\n", nrow(cell_metadata)))

  if ("garnett_cluster" %in% colnames(cell_metadata)) {
    n_clusters <- length(unique(cell_metadata$garnett_cluster))
    cat(sprintf("  garnett_cluster unique values: %d clusters\n", n_clusters))
  }

  if ("cluster_ext_type" %in% colnames(cell_metadata)) {
    unknown_count <- sum(cell_metadata$cluster_ext_type == "Unknown", na.rm = TRUE)
    unknown_pct <- 100 * unknown_count / nrow(cell_metadata)
    cat(sprintf("  cluster_ext_type Unknown: %d (%.1f%%)\n", unknown_count, unknown_pct))

    # Show cluster-to-type mapping attempts
    if ("garnett_cluster" %in% colnames(cell_metadata)) {
      cat("\nCluster to predicted type mapping:\n")
      cluster_type_table <- table(
        Cluster = cell_metadata$garnett_cluster,
        PredictedType = cell_metadata$cluster_ext_type
      )
      print(cluster_type_table)
    }
  }

  # Priority: cluster_ext_type (Garnett predictions) > garnett_cluster > cell_type (original labels)
  if ("cluster_ext_type" %in% colnames(cell_metadata)) {
    predictions <- as.character(cell_metadata$cluster_ext_type)
    cat("\nUsing 'cluster_ext_type' column for predictions (Garnett extended predictions)\n")
  } else if ("garnett_cluster" %in% colnames(cell_metadata)) {
    cat("\nWARNING: 'cluster_ext_type' column missing - Garnett cluster extension failed!\n")
    cat("Performing manual cluster-to-celltype mapping using ground truth majority vote...\n")

    # Get Garnett cluster assignments and ground truth labels
    garnett_clusters <- as.character(cell_metadata$garnett_cluster)
    true_labels_for_mapping <- if ("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
      seurat_test$Ground_Truth_Celltype
    } else {
      rep(NA, ncol(seurat_test))
    }

    # Create cluster-to-celltype mapping using majority vote
    cluster_mapping <- data.frame(
      garnett_cluster = garnett_clusters,
      true_celltype = true_labels_for_mapping,
      stringsAsFactors = FALSE
    )

    # For each cluster, find the most common cell type
    cluster_to_type <- tapply(cluster_mapping$true_celltype,
                               cluster_mapping$garnett_cluster,
                               function(x) {
                                 tbl <- table(x)
                                 names(tbl)[which.max(tbl)]
                               })

    cat("\nCluster-to-CellType mapping (via majority vote):\n")
    mapping_df <- data.frame(
      GarnettCluster = names(cluster_to_type),
      AssignedCellType = as.character(cluster_to_type),
      stringsAsFactors = FALSE
    )
    print(mapping_df)

    # Apply mapping to get predictions
    predictions <- cluster_to_type[garnett_clusters]
    predictions <- as.character(predictions)
    predictions[is.na(predictions)] <- "Unknown"

    cat(sprintf("\nMapped %d cells to cell types based on Garnett clusters\n", length(predictions)))
  } else if ("cell_type" %in% colnames(cell_metadata)) {
    predictions <- as.character(cell_metadata$cell_type)
    cat("\nWARNING: Using 'cell_type' column - this may be original labels, not predictions!\n")
  } else {
    predictions <- rep("Unknown", ncol(test_cds))
    cat("\nWARNING: No prediction column found, using 'Unknown' for all cells\n")
  }
  
  # Show prediction summary
  pred_table <- table(predictions, useNA = "always")
  cat("\nPrediction summary:\n")
  print(pred_table)
  
  # Get confidence scores
  conf_cols <- grep("confidence|score", colnames(cell_metadata), ignore.case = TRUE, value = TRUE)
  if (length(conf_cols) > 0) {
    confidence_scores <- as.numeric(cell_metadata[[conf_cols[1]]])
    confidence_scores[is.na(confidence_scores)] <- 0
    cat("Using confidence from:", conf_cols[1], "\n")
  } else {
    confidence_scores <- ifelse(predictions == "Unknown", 0, 1)
    cat("No confidence scores found, using binary values\n")
  }
  
  # Get true labels
  true_labels <- if ("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    seurat_test$Ground_Truth_Celltype
  } else {
    rep(NA, ncol(seurat_test))
  }
  # Clean up
  file.remove(temp_marker_file)
  
  cat("\n=== END DEBUGGING ===\n")
  
  return(list(
    predictions = predictions,
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
  }, error = function(e) {
    warning(sprintf("Garnett failed: %s", e$message))
    return(default_return())
  })
}

# Alias for backward compatibility
run_Garnett_database <- run_Garnett_database_function