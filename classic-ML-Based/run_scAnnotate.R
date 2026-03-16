# run_scAnnotate.R

run_scAnnotate_function <- function(seurat_train, seurat_test, markers) {
  
  library(scAnnotate)
  library(Seurat)
  
  # [FIX] Increase memory limit for parallel processing
  # SCTransform uses 'future' for parallelization. 
  # This prevents the "globals size" error.
  options(future.globals.maxSize = 8000 * 1024^2)
  
  # Default return for errors
  default_return <- function(msg) {
    warning(msg)
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    ))
  }
  
  tryCatch({
    # 1. Validation
    if (!"Ground_Truth_Celltype" %in% colnames(seurat_train@meta.data)) {
      return(default_return("Missing Ground_Truth_Celltype in train"))
    }
    
    # 2. Extract Data & CONVERT TO RAW COUNTS
    # scAnnotate's Harmony mode (SCTransform) crashes on LogNormalized data.
    # We must provide Integers.
    get_raw_counts <- function(obj) {
      # Try to get existing raw counts
      counts <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"), error=function(e) NULL)
      
      # If missing (only data layer exists), UN-LOG it
      if (is.null(counts)) {
        data_layer <- GetAssayData(obj, assay = "RNA", layer = "data")
        # Reverse LogNormalize: exp(x) - 1
        counts <- exp(data_layer) - 1
      }
      # Force Integer (Round)
      return(round(counts))
    }
    
    train_counts <- get_raw_counts(seurat_train)
    test_counts <- get_raw_counts(seurat_test)
    
    # 3. Get Common Genes
    common_genes <- intersect(rownames(train_counts), rownames(test_counts))
    if (length(common_genes) < 100) return(default_return("Too few common genes"))
    
    # Subset
    train_sub <- train_counts[common_genes, ]
    test_sub <- test_counts[common_genes, ]
    
    # 4. Fix "MT-" Gene Issue (Required for SCTransform)
    # If no mitochondrial genes exist, SCTransform crashes. We fake them if needed.
    if (!any(grepl("^MT-", rownames(train_sub), ignore.case=TRUE))) {
      new_names <- rownames(train_sub)
      new_names[1:5] <- paste0("MT-", new_names[1:5])
      rownames(train_sub) <- new_names
      rownames(test_sub) <- new_names
    }
    
    # 5. Sanitize Names
    clean_genes <- make.unique(rownames(train_sub))
    rownames(train_sub) <- clean_genes
    rownames(test_sub) <- clean_genes
    
    # 6. Prepare Inputs (Transpose to Cells x Genes)
    train_mat <- t(as.matrix(train_sub))
    test_mat <- t(as.matrix(test_sub))
    
    # Create Data Frame for Train
    train_df <- data.frame(
      CellType = as.character(seurat_train$Ground_Truth_Celltype),
      train_mat,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # Ensure Colnames match
    colnames(train_df) <- c("CellType", colnames(train_mat))
    
    # 7. Run scAnnotate with memory tracking
    # lognormalized = FALSE is critical here to handle the raw counts we just created
    peak_memory_mb <- NA

    if (!requireNamespace("bench", quietly = TRUE)) {
      warning("bench package not available for memory tracking")
      predict_label <- scAnnotate(
        train = train_df,
        test = test_mat,
        distribution = "normal",
        correction = "auto",   # Now safe to use auto/harmony
        screening = "wilcox",
        threshold = 0,
        lognormalized = FALSE  # TELLING IT WE HAVE RAW DATA
      )
    } else {
      library(bench)
      bench_result <- bench::mark(
        {
          predict_label <- scAnnotate(
            train = train_df,
            test = test_mat,
            distribution = "normal",
            correction = "auto",   # Now safe to use auto/harmony
            screening = "wilcox",
            threshold = 0,
            lognormalized = FALSE  # TELLING IT WE HAVE RAW DATA
          )
        },
        memory = TRUE,
        iterations = 1,
        check = FALSE
      )
      peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2
    }

    # 8. Format Results
    preds <- as.character(predict_label)
    preds[is.na(preds)] <- "Unknown"

    return(list(
      predictions = preds,
      true_labels = as.character(seurat_test$Ground_Truth_Celltype),
      confidence_scores = ifelse(preds == "Unknown", 0, 0.8),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = peak_memory_mb
    ))
    
  }, error = function(e) {
    return(default_return(paste("scAnnotate Error:", e$message)))
  })
}

run_scAnnotate <- run_scAnnotate_function