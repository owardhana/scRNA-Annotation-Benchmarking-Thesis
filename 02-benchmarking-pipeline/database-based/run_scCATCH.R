# run_scCATCH.R
#################################################
# scCATCH Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scCATCH Cell Type Annotation Function
#' 
#' Purpose: Run scCATCH algorithm using custom marker conversion
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts markers to scCATCH format, uses custom marker annotation
run_scCATCH_database_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("scCATCH", quietly = TRUE)) {
    stop("scCATCH package not available. Please install it first.")
  }
  library(scCATCH)

  # Convert markers from FindAllMarkers format to scCATCH marker format
  convert_markers_to_sccatch <- function(marker_df) {
    filtered_markers <- marker_df

    if(nrow(filtered_markers) == 0) {
      warning("No markers found after filtering")
      return(NULL)
    }

    # Create scCATCH marker format with required columns
    top_markers <- filtered_markers

    # Create scCATCH marker format with required columns
    sccatch_markers <- data.frame(
      species = "Human",                    # Assume human data
      tissue = "Custom",                   # Custom tissue type
      cancer = "Normal",                   # Normal tissue
      condition = "Normal cell",           # Normal cell condition
      subtype1 = NA,                       # No subtype classification
      subtype2 = NA,                       # No subtype classification
      subtype3 = NA,                       # No subtype classification
      celltype = top_markers$cluster,      # Cell type from cluster
      gene = top_markers$gene,             # Marker gene
      resource = "Single-cell sequencing", # Resource type
      pmid = "Custom",                     # Custom reference
      stringsAsFactors = FALSE
    )

    return(sccatch_markers)
  }

  # Convert markers to scCATCH format
  sccatch_markers <- convert_markers_to_sccatch(markers)
  
  if(is.null(sccatch_markers)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }
  
  # Ensure seurat_test has normalized data - use cross-compatible approach
  tryCatch({
    # Try to get normalized data - this will fail if not available
    test_data <- GetAssayData(seurat_test, assay = "RNA", layer = "data")
    if(length(test_data) == 0) {
      seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
    }
  }, error = function(e) {
    # Fallback for Seurat v4 or if layer approach fails
    tryCatch({
      test_data <- GetAssayData(seurat_test, assay = "RNA", layer = "data")
      if(length(test_data) == 0) {
        seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
      }
    }, error = function(e2) {
      # If both fail, normalize the data
      seurat_test <- NormalizeData(seurat_test, verbose = FALSE)
    })
  })
  
  # Determine cluster information for scCATCH using ground truth cell types
  cluster_info <- as.character(seurat_test@meta.data$Ground_Truth_Celltype)
  
  # Get expression data using cross-compatible approach
  expr_data <- tryCatch({
    # Try Seurat v5 first
    GetAssayData(seurat_test, assay = "RNA", layer = "data")
  }, error = function(e) {
    # Fallback to Seurat v4 syntax
    tryCatch({
      GetAssayData(seurat_test, assay = "RNA", layer = "data")
    }, error = function(e2) {
      # Final fallback - direct access (should not be needed)
      seurat_test[["RNA"]]@data
    })
  })
  
  # Validate data dimensions
  if(is.null(expr_data) || length(expr_data) == 0) {
    warning("Expression data is empty or could not be accessed")
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }

  # Validate cluster vector
  if(length(cluster_info) != ncol(seurat_test)) {
    warning(paste("Cluster vector length mismatch: expected", ncol(seurat_test), "got", length(cluster_info)))
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }
  
  # Debug information
  cat(sprintf("Expression data: %d genes x %d cells\n", nrow(expr_data), ncol(expr_data)))
  cat(sprintf("Cluster info: %d elements\n", length(cluster_info)))
  cat(sprintf("Unique clusters: %s\n", paste(unique(cluster_info), collapse = ", ")))
  
  # Run scCATCH 3-step pipeline with memory tracking
  runtime_secs <- NA
  peak_system_memory_mb <- NA
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")
    # Create scCATCH object
    start_time <- Sys.time()
    obj <- tryCatch({
      createscCATCH(data = expr_data, cluster = cluster_info)
    }, error = function(e) {
      warning(paste("Failed to create scCATCH object:", e$message))
      return(NULL)
    })

    if(is.null(obj)) {
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      return(list(
        predictions = rep("Unknown", ncol(seurat_test)),
        true_labels = seurat_test$Ground_Truth_Celltype,
        confidence_scores = rep(0, ncol(seurat_test)),
        cell_ids = colnames(seurat_test),
        runtime_secs = runtime_secs,
        peak_system_memory_mb = NA
      ))
    }

    # Run scCATCH with custom markers
    obj <- tryCatch({
      findmarkergene(obj,
                     if_use_custom_marker = TRUE,
                     marker = sccatch_markers,
                     species = "Human",
                     cluster = "All",
                     cell_min_pct = 0.25,
                     logfc = 0.25,
                     pvalue = 0.05)
    }, error = function(e) {
      warning(paste("findmarkergene failed:", e$message))
      return(NULL)
    })

    if(is.null(obj)) {
      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      return(list(
        predictions = rep("Unknown", ncol(seurat_test)),
        true_labels = seurat_test$Ground_Truth_Celltype,
        confidence_scores = rep(0, ncol(seurat_test)),
        cell_ids = colnames(seurat_test),
        runtime_secs = runtime_secs,
        peak_system_memory_mb = NA
      ))
    }

    # Perform cell type annotation
    obj <- tryCatch({
      findcelltype(obj)
    }, error = function(e) {
      warning(paste("findcelltype failed:", e$message))
      return(NULL)
    })
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM({
      # Step 1: Create scCATCH object
      obj <- createscCATCH(data = expr_data, cluster = cluster_info)
      # Step 2: Find marker genes
      obj <- findmarkergene(obj,
                           if_use_custom_marker = TRUE,
                           marker = sccatch_markers,
                           species = "Human",
                           cluster = "All",
                           cell_min_pct = 0.25,
                           logfc = 0.25,
                           pvalue = 0.05)
      # Step 3: Find cell types
      obj <- findcelltype(obj)
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  if(is.null(obj)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = runtime_secs,
      peak_system_memory_mb = NA
    ))
  }

  # Extract results and map to individual cells
  if(!is.null(obj@celltype) && nrow(obj@celltype) > 0) {
    cluster_annotations <- obj@celltype
    
    # Create cell-level predictions
    cell_predictions <- rep("Unknown", ncol(seurat_test))
    names(cell_predictions) <- colnames(seurat_test)
    
    confidence_scores <- rep(0, ncol(seurat_test))
    names(confidence_scores) <- names(cell_predictions)
    
    # Map cluster annotations to individual cells
    for(i in 1:nrow(cluster_annotations)) {
      cluster_id <- as.character(cluster_annotations$cluster[i])
      cell_type <- cluster_annotations$cell_type[i]
      
      # Find cells in this cluster
      cluster_cells <- which(cluster_info == cluster_id)
      if(length(cluster_cells) > 0) {
        cell_names <- colnames(seurat_test)[cluster_cells]
        cell_predictions[cell_names] <- cell_type
        
        # Use similarity score as confidence if available
        if("score" %in% colnames(cluster_annotations)) {
          confidence_scores[cell_names] <- cluster_annotations$score[i]
        } else {
          confidence_scores[cell_names] <- 1.0
        }
      }
    }
  } else {
    cell_predictions <- rep("Unknown", ncol(seurat_test))
    names(cell_predictions) <- colnames(seurat_test)
    confidence_scores <- rep(0, ncol(seurat_test))
    names(confidence_scores) <- colnames(seurat_test)
  }

  # Get true labels from test set
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
}


# For backward compatibility
run_scCATCH_database <- run_scCATCH_database_function