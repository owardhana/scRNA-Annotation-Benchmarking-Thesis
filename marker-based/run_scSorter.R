# run_scSorter.R
#################################################
# scSorter Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scSorter Cell Type Annotation Function
#' 
#' Purpose: Run scSorter algorithm on test data using marker signatures
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict (already normalized with variable features)
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Minimal preprocessing, convert markers to scSorter format, run scSorter
run_scSorter_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("scSorter", quietly = TRUE)) {
    stop("scSorter package not available. Please install it first.")
  }
  library(scSorter)

  # Convert markers from FindAllMarkers format to scSorter annotation format
  convert_markers_to_scsorter <- function(marker_df, markers_per_type = 5) { #instead of 50 
    # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
    filtered_markers <- marker_df[marker_df$avg_log2FC >= 1 & #instead of 0.5
                                 marker_df$p_val_adj < 0.25 & #instead of 0.05
                                 marker_df$pct.1 >= 0.25, ] #instead of 0.15

    # Take top N markers per cell type (standardized: 50)
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr required for marker selection")
    }
    library(dplyr)

    filtered_markers <- filtered_markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::slice_head(n = markers_per_type) %>%
      dplyr::ungroup()

    # scSorter expects: Type, Marker, Weight columns
    anno <- data.frame(
      Type = filtered_markers$cluster,
      Marker = filtered_markers$gene,
      Weight = filtered_markers$avg_log2FC,
      stringsAsFactors = FALSE
    )
    return(anno)
  }

  # Simplified preprocessing since seurat_test is already normalized with variable features
  # Get pre-computed variable features
  topgenes <- VariableFeatures(seurat_test)
  if(length(topgenes) == 0) {
    warning("No variable features found, using all genes")
    topgenes <- rownames(seurat_test)
  }
  
  # Get normalized expression data
  expr <- GetAssayData(seurat_test, layer = "data")
  
  # Filter for highly variable genes expressed in >10% of cells
  topgene_filter <- rowSums(as.matrix(expr)[topgenes, ] != 0) > ncol(expr) * 0.1
  topgenes <- topgenes[topgene_filter]
  
  # Convert markers to scSorter annotation format
  anno <- convert_markers_to_scsorter(markers)

  # Filter annotation to only include genes present in expression data
  anno_filtered <- anno[anno$Marker %in% rownames(expr), ]
  
  if(nrow(anno_filtered) == 0) {
    warning("No marker genes found in expression data")
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test)
    ))
  }
  
  # Include marker genes and top variable genes, then subset expression data
  picked_genes <- unique(c(anno_filtered$Marker, topgenes))
  expr <- expr[rownames(expr) %in% picked_genes, ]

  # Run scSorter with memory tracking
  peak_memory_mb <- NA
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("bench package not available for memory tracking")
    results <- tryCatch({
      scSorter(expr, anno_filtered)
    }, error = function(e) {
      warning(paste("scSorter failed:", e$message))
      return(NULL)
    })
  } else {
    library(bench)
    bench_result <- bench::mark(
      {
        results <- scSorter(expr, anno_filtered)
      },
      memory = TRUE,
      iterations = 1,
      check = FALSE
    )
    # Extract peak memory in MB
    peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
  }

  if(is.null(results)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }

  # Extract results
  cell_predictions <- results$Pred_Type
  cell_scores <- results$Scores
  
  # Get confidence scores (maximum score for each cell)
  if(!is.null(cell_scores)) {
    confidence_scores <- apply(cell_scores, 1, max, na.rm = TRUE)
  } else {
    # If no scores available, set confidence to 1 for all predictions
    confidence_scores <- rep(1, length(cell_predictions))
  }
  
  # Get true labels from test set
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    peak_memory_mb = peak_memory_mb
  ))
}

# For backward compatibility
run_scSorter <- run_scSorter_function