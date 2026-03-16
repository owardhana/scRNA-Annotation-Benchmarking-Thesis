# run_scType.R
#################################################
# scType Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scType Cell Type Annotation Function
#' 
#' Purpose: Run scType algorithm using official wrapper with custom marker conversion
#' Inputs:
#'   - seurat_train: Training Seurat object (not used directly, for interface consistency)
#'   - seurat_test: Test Seurat object to predict (will be scaled if needed)
#'   - markers: Marker genes dataframe from FindAllMarkers()
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Converts markers to scType Excel format, uses official run_sctype wrapper
run_scType_function <- function(seurat_train, seurat_test, markers) {
  
  # Load all required libraries for scType
  required_packages <- c("dplyr", "Seurat", "HGNChelper", "openxlsx")
  
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "required for scType. Please install it with: install.packages('", pkg, "')"))
    }
    library(pkg, character.only = TRUE)
  }
  
  # Load all required scType source files in correct order
  tryCatch({
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
  }, error = function(e) {
    stop(paste("Failed to load scType source files:", e$message))
  })
  
  # Verify key functions are loaded
  required_functions <- c("gene_sets_prepare", "sctype_score", "run_sctype")
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  if(length(missing_functions) > 0) {
    stop(paste("Missing scType functions:", paste(missing_functions, collapse = ", ")))
  }

  # Convert FindAllMarkers output to scType Excel format
  convert_to_sctype_excel <- function(markers) {
    
    # Don't pre-process gene symbols - let scType handle them natively
    # This avoids duplicate processing and warnings
    
    cat("Converting", nrow(markers), "markers to scType format...\n")

    # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
    filtered_markers <- markers[markers$avg_log2FC >= 0.5 &
                               markers$p_val_adj < 0.05 &
                               markers$pct.1 >= 0.15, ]

    if(nrow(filtered_markers) == 0) {
      warning("No high-quality markers found")
      return(NULL)
    }

    # Take top 50 markers per cell type (standardized)
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr required for marker selection")
    }
    library(dplyr)

    filtered_markers <- filtered_markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::ungroup()

    # Group by cell type and prepare scType format
    cell_types <- unique(filtered_markers$cluster)
    cat("Found", length(cell_types), "cell types:\n")
    
    # Try to use the template file first for proper format
    template_file <- "misc/ScTypeDB_full.xlsx"
    sctype_data <- NULL
    
    if(file.exists(template_file)) {
      tryCatch({
        template_data <- read.xlsx(template_file)
        cat("Using template file with", nrow(template_data), "rows and columns:", paste(colnames(template_data), collapse=", "), "\n")
        
        # Create custom data with same structure as template
        sctype_data <- data.frame(
          tissueType = character(),
          cellName = character(),
          geneSymbolmore1 = character(),
          geneSymbolmore2 = character(),
          shortName = character(),
          stringsAsFactors = FALSE
        )
        
        # Template has exactly 5 columns, so structure is already correct
        
      }, error = function(e) {
        cat("Failed to use template:", e$message, "\n")
        sctype_data <<- data.frame()
      })
    } else {
      cat("Template file not found, using manual format\n")
      sctype_data <- data.frame(
        tissueType = character(),
        cellName = character(),
        geneSymbolmore1 = character(),
        geneSymbolmore2 = character(),
        shortName = character(),
        stringsAsFactors = FALSE
      )
    }
    
    # Build the data
    for(ct in cell_types) {
      ct_markers <- filtered_markers[filtered_markers$cluster == ct, ]
      
      # Get positive markers (genes upregulated in this cell type)
      pos_markers <- ct_markers$gene
      
      # Only include if we have positive markers
      if(length(pos_markers) > 0) {
        # Get negative markers (strong markers from other cell types) - optional
        neg_markers <- filtered_markers$gene[filtered_markers$cluster != ct &
                                             filtered_markers$avg_log2FC >= 2.0]
        
        # Create scType format row with exact template structure
        new_row <- data.frame(
          tissueType = "Custom",  # Use Custom tissue type to match known_tissue_type parameter
          cellName = as.character(ct),
          geneSymbolmore1 = paste(pos_markers, collapse = ","),
          geneSymbolmore2 = if(length(neg_markers) > 0) paste(head(neg_markers, 10), collapse = ",") else "",
          shortName = as.character(ct),  # Add shortName column matching template
          stringsAsFactors = FALSE
        )
        
        # Template structure is exactly 5 columns, so new_row already matches
        
        sctype_data <- rbind(sctype_data, new_row)
        
        cat("  ", ct, ":", length(pos_markers), "positive,", length(neg_markers), "negative markers\n")
      }
    }
    
    # Validate the data before returning
    if(nrow(sctype_data) == 0) {
      warning("No valid cell types with markers found")
      return(NULL)
    }
    
    # Ensure no empty gene lists
    sctype_data <- sctype_data[sctype_data$geneSymbolmore1 != "", ]
    
    if(nrow(sctype_data) == 0) {
      warning("All cell types have empty marker lists")
      return(NULL)
    }
    
    return(sctype_data)
  }

  # Convert markers to scType Excel format
  sctype_data <- convert_to_sctype_excel(markers)
  
  if(is.null(sctype_data)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test)
    ))
  }
  
  # Debug the final data structure
  cat("\nFinal scType data structure:\n")
  str(sctype_data)
  cat("\nColumn names:", paste(colnames(sctype_data), collapse=", "), "\n")
  cat("Sample rows:\n")
  if(nrow(sctype_data) > 0) {
    print(head(sctype_data, 2))
  }
  
  # Create temporary Excel file for scType
  temp_file <- tempfile(pattern = "sctype_markers_", fileext = ".xlsx")
  tryCatch({
    # Write Excel file with correct format - use same sheet name as template
    write.xlsx(sctype_data, temp_file, sheetName = "Sheet1", rowNames = FALSE)
    
    # Validate file creation
    if(!file.exists(temp_file) || file.size(temp_file) == 0) {
      stop("Failed to create valid Excel file")
    }
    
    cat("\nCreated marker file:", temp_file, "(size:", file.size(temp_file), "bytes)\n")
    cat("Cell types in file:", nrow(sctype_data), "\n")
    
    # Validate by reading back
    test_read <- read.xlsx(temp_file)
    cat("Read-back validation: ", nrow(test_read), "rows x", ncol(test_read), "columns\n")
    cat("Read-back columns:", paste(colnames(test_read), collapse=", "), "\n")
    
  }, error = function(e) {
    stop(paste("Failed to create marker file:", e$message))
  })
  
  # Ensure seurat_test has scaled data
  if(!"scale.data" %in% names(seurat_test@assays$RNA@layers) || 
     length(seurat_test@assays$RNA@layers$scale.data) == 0) {
    seurat_test <- ScaleData(seurat_test, features = rownames(seurat_test), verbose = FALSE)
  }

  # Run scType using official wrapper with memory tracking
  peak_memory_mb <- NA

  # Ensure we have clusters for scType (it expects clustered data)
  if(!"seurat_clusters" %in% colnames(seurat_test@meta.data)) {
    # Run basic clustering if not present
    cat("Adding clustering to seurat_test...\n")
    seurat_test <- FindNeighbors(seurat_test, verbose = FALSE)
    seurat_test <- FindClusters(seurat_test, verbose = FALSE)
  }

  results <- tryCatch({
    cat("Running scType with file:", temp_file, "\n")
    cat("Seurat object:", ncol(seurat_test), "cells x", nrow(seurat_test), "genes\n")

    if (!requireNamespace("bench", quietly = TRUE)) {
      warning("bench package not available for memory tracking")
      # Suppress ALL warnings from scType and HGNChelper
      suppressMessages(suppressWarnings({
        run_sctype(seurat_test,
                   assay = "RNA",
                   scaled = TRUE,
                   known_tissue_type = "Custom",
                   custom_marker_file = temp_file,
                   name = "sctype_classification")
      }))
    } else {
      library(bench)
      bench_result <- bench::mark(
        {
          result <- suppressMessages(suppressWarnings({
            run_sctype(seurat_test,
                       assay = "RNA",
                       scaled = TRUE,
                       known_tissue_type = "Custom",
                       custom_marker_file = temp_file,
                       name = "sctype_classification")
          }))
        },
        memory = TRUE,
        iterations = 1,
        check = FALSE
      )
      # Extract peak memory in MB
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
      result
    }
  }, error = function(e) {
    cat("scType wrapper failed:", e$message, "\n")
    
    # Try fallback approach using direct sctype_score
    cat("Attempting fallback approach with direct sctype_score...\n")
    fallback_result <- tryCatch({

      # Convert markers to gene sets manually
      # STANDARDIZED FILTERING: avg_log2FC >= 0.5, p_val_adj < 0.05, pct.1 >= 0.15
      fallback_markers <- markers[markers$avg_log2FC >= 0.5 &
                                 markers$p_val_adj < 0.05 &
                                 markers$pct.1 >= 0.15, ]

      # Take top 50 markers per cell type (standardized)
      fallback_markers <- fallback_markers %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        dplyr::slice_head(n = 50) %>%
        dplyr::ungroup()

      gs_positive <- list()
      for(ct in unique(fallback_markers$cluster)) {
        ct_markers <- fallback_markers$gene[fallback_markers$cluster == ct]
        if(length(ct_markers) > 0) {
          gs_positive[[ct]] <- as.character(ct_markers)
        }
      }

      cat("Created gene sets for", length(gs_positive), "cell types\n")

      # Use sctype_score directly with memory tracking
      # gene_names_to_uppercase = FALSE: avoid case mismatch with custom/synthetic gene names
      if (!requireNamespace("bench", quietly = TRUE)) {
        es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(seurat_test, assay="RNA", layer="scale.data")),
                               scaled = TRUE,
                               gs = gs_positive,
                               gene_names_to_uppercase = FALSE)
      } else {
        library(bench)
        bench_result_fallback <- bench::mark(
          {
            es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(seurat_test, assay="RNA", layer="scale.data")),
                                   scaled = TRUE,
                                   gs = gs_positive,
                                   gene_names_to_uppercase = FALSE)
          },
          memory = TRUE,
          iterations = 1,
          check = FALSE
        )
        # Extract peak memory in MB
        peak_memory_mb <<- as.numeric(bench_result_fallback$mem_alloc) / 1024^2
      }

      # Extract predictions from sctype_score results
      # es.max: rows = cell types, columns = cells
      cell_predictions <- rownames(es.max)[apply(es.max, 2, which.max)]
      confidence_scores <- apply(es.max, 2, max)

      # Create mock results object with proper cell ID alignment
      mock_results <- seurat_test  # Keep original object unchanged

      # Align results with original Seurat object cell IDs using match
      mock_results@meta.data$sctype_classification <- cell_predictions[match(colnames(seurat_test), colnames(es.max))]
      mock_results@meta.data$sctype_scores <- confidence_scores[match(colnames(seurat_test), colnames(es.max))]

      cat("Fallback approach succeeded\n")
      return(mock_results)

    }, error = function(e2) {
      cat("Fallback approach also failed:", e2$message, "\n")
      return(NULL)
    })
    
    if(!is.null(fallback_result)) {
      return(fallback_result)
    } else {
      warning("Both scType wrapper and fallback approach failed")
      return(NULL)
    }
    
  }, finally = {
    # Clean up temporary file
    if(file.exists(temp_file)) {
      file.remove(temp_file)
    }
  })
  
  if(is.null(results)) {
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }

  # Extract predictions from scType results
  # scType adds the classification to metadata
  if("sctype_classification" %in% colnames(results@meta.data)) {
    cell_predictions <- results@meta.data$sctype_classification
  } else {
    # Try alternative column names
    sctype_cols <- grep("sctype", colnames(results@meta.data), ignore.case = TRUE, value = TRUE)
    if(length(sctype_cols) > 0) {
      cell_predictions <- results@meta.data[[sctype_cols[1]]]
      cat("Using scType column:", sctype_cols[1], "\n")
    } else {
      # List available columns for debugging
      cat("Available metadata columns:", paste(colnames(results@meta.data), collapse = ", "), "\n")
      warning("Could not find scType predictions in results")
      cell_predictions <- rep("Unknown", ncol(results))
    }
  }
  
  # Debug prediction summary
  cat("scType predictions summary:\n")
  print(table(cell_predictions, useNA = "ifany"))
  
  # Get confidence scores from scType scores if available
  if("sctype_scores" %in% colnames(results@meta.data)) {
    confidence_scores <- results@meta.data$sctype_scores
  } else {
    # Try to find score columns
    score_cols <- grep("score", colnames(results@meta.data), ignore.case = TRUE, value = TRUE)
    if(length(score_cols) > 0) {
      confidence_scores <- results@meta.data[[score_cols[1]]]
    } else {
      # If no confidence scores, set to 1 for all predictions except "Unknown"
      confidence_scores <- ifelse(cell_predictions == "Unknown", 0, 1)
    }
  }
  
  # Handle unknown/unassigned predictions
  cell_predictions[is.na(cell_predictions) | cell_predictions == "" | cell_predictions == "Unknown"] <- "Unknown"
  
  # Ensure confidence scores are numeric and same length as predictions
  confidence_scores <- as.numeric(confidence_scores)
  confidence_scores[is.na(confidence_scores)] <- 0
  
  # Ensure lengths match
  if(length(confidence_scores) != length(cell_predictions)) {
    warning("Confidence scores length mismatch, using default values")
    confidence_scores <- ifelse(cell_predictions == "Unknown", 0, 1)
  }
  
  # Final validation
  if(length(cell_predictions) != ncol(seurat_test)) {
    warning("Prediction length mismatch with input cells")
    # Pad or truncate as needed
    if(length(cell_predictions) < ncol(seurat_test)) {
      cell_predictions <- c(cell_predictions, rep("Unknown", ncol(seurat_test) - length(cell_predictions)))
      confidence_scores <- c(confidence_scores, rep(0, ncol(seurat_test) - length(confidence_scores)))
    } else {
      cell_predictions <- cell_predictions[1:ncol(seurat_test)]
      confidence_scores <- confidence_scores[1:ncol(seurat_test)]
    }
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
run_scType <- run_scType_function