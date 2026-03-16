# run_SingleCellNet.R
#################################################
# SingleCellNet Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' SingleCellNet Cell Type Annotation Function
#' 
#' Purpose: Run SingleCellNet algorithm using random forest-based classification
#' Inputs:
#'   - seurat_train: Training Seurat object (used to train classifier)
#'   - seurat_test: Test Seurat object to predict 
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses SingleCellNet's scn_train and scn_predict workflow
run_singleCellNet_function <- function(seurat_train, seurat_test, markers) {
  
  # Load required library
  if (!requireNamespace("singleCellNet", quietly = TRUE)) {
    stop("singleCellNet package not available. Please install it first.")
  }
  library(singleCellNet)
  library(Seurat)

  # Default return with memory field
  default_return <- function() {
    list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = seurat_test$Ground_Truth_Celltype,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      peak_memory_mb = NA
    )
  }

  tryCatch({

  # Extract data from Seurat objects manually (extractSeurat uses Seurat v4 slot API)
  train_data <- tryCatch({
    list(
      sampTab = seurat_train@meta.data,
      expDat = as.matrix(GetAssayData(seurat_train, assay = "RNA", layer = "counts"))
    )
  }, error = function(e) {
    warning(paste("Failed to extract training data:", e$message))
    return(NULL)
  })

  test_data <- tryCatch({
    list(
      sampTab = seurat_test@meta.data,
      expDat = as.matrix(GetAssayData(seurat_test, assay = "RNA", layer = "counts"))
    )
  }, error = function(e) {
    warning(paste("Failed to extract test data:", e$message))
    return(NULL)
  })

  # Ensure proper column names in sample tables
  if(!is.null(train_data) && !is.null(train_data$sampTab)) {
    if(!"Ground_Truth_Celltype" %in% colnames(train_data$sampTab)) {
      train_data$sampTab$Ground_Truth_Celltype <- seurat_train$Ground_Truth_Celltype
    }
    if(!"cell_id" %in% colnames(train_data$sampTab)) {
      train_data$sampTab$cell_id <- rownames(train_data$sampTab)
    }
  }

  if(!is.null(test_data) && !is.null(test_data$sampTab)) {
    if(!"Ground_Truth_Celltype" %in% colnames(test_data$sampTab)) {
      test_data$sampTab$Ground_Truth_Celltype <- seurat_test$Ground_Truth_Celltype
    }
    if(!"cell_id" %in% colnames(test_data$sampTab)) {
      test_data$sampTab$cell_id <- rownames(test_data$sampTab)
    }
  }

  if(is.null(train_data) || is.null(test_data)) {
    return(default_return())
  }

  # Get common genes between training and test data
  commonGenes <- intersect(rownames(train_data$expDat), rownames(test_data$expDat))

  if(length(commonGenes) < 100) {
    warning("Too few common genes between training and test data")
    return(default_return())
  }

  expTrain <- train_data$expDat[commonGenes, ]
  expTest <- test_data$expDat[commonGenes, ]

  # Track peak memory usage
  peak_memory_mb <- NA

  # Run singleCellNet training and prediction with memory tracking
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("bench package not available for memory tracking")

    # Train classifier
    set.seed(100)
    class_info <- tryCatch({
      scn_train(stTrain = train_data$sampTab,
                expTrain = expTrain,
                dLevel = "Ground_Truth_Celltype",
                colName_samp = "cell_id",
                nTopGenes = 10,
                nRand = 70,
                nTrees = 1000,
                nTopGenePairs = 25)
    }, error = function(e) {
      warning(paste("Failed to train SingleCellNet classifier:", e$message))
      return(NULL)
    })

    if(is.null(class_info)) {
      return(default_return())
    }

    # Predict on test data
    classRes_test <- tryCatch({
      scn_predict(class_info[['cnProc']], expTest, nrand = 70)
    }, error = function(e) {
      warning(paste("Failed to predict with SingleCellNet:", e$message))
      return(NULL)
    })

    if(is.null(classRes_test)) {
      return(default_return())
    }

    # Get predicted cell types using get_cate function
    stTest_pred <- tryCatch({
      get_cate(classRes_test, test_data$sampTab,
               dLevel = "Ground_Truth_Celltype",
               sid = "cell_id",
               nrand = 70)
    }, error = function(e) {
      warning(paste("Failed to extract categories:", e$message))
      return(NULL)
    })

    if(is.null(stTest_pred)) {
      return(default_return())
    }
  } else {
    library(bench)
    bench_result <- bench::mark(
      {
        # Train classifier
        set.seed(100)
        class_info <- scn_train(stTrain = train_data$sampTab,
                                expTrain = expTrain,
                                dLevel = "Ground_Truth_Celltype",
                                colName_samp = "cell_id",
                                nTopGenes = 10,
                                nRand = 70,
                                nTrees = 1000,
                                nTopGenePairs = 25)

        # Predict on test data
        classRes_test <- scn_predict(class_info[['cnProc']], expTest, nrand = 70)

        # Get predicted cell types using get_cate function
        stTest_pred <- get_cate(classRes_test, test_data$sampTab,
                                dLevel = "Ground_Truth_Celltype",
                                sid = "cell_id",
                                nrand = 70)
      },
      memory = TRUE,
      iterations = 1,
      check = FALSE
    )
    peak_memory_mb <- as.numeric(bench_result$mem_alloc) / 1024^2

    # Check for errors in bench block
    if(!exists("stTest_pred") || is.null(stTest_pred)) {
      warning("singleCellNet training or prediction failed")
      return_val <- default_return()
      return_val$peak_memory_mb <- peak_memory_mb
      return(return_val)
    }
  }

  # Extract predictions and true labels
  cell_predictions <- stTest_pred$category
  true_labels <- seurat_test$Ground_Truth_Celltype
  
  # Extract confidence scores from classification results
  # SingleCellNet provides classification scores matrix (cells x cell_types)
  if(!is.null(classRes_test) && is.matrix(classRes_test)) {
    # Get maximum probability for each cell as confidence
    confidence_scores <- apply(classRes_test, 1, max, na.rm = TRUE)
    # Normalize scores to 0-1 range if needed
    if(max(confidence_scores, na.rm = TRUE) > 1) {
      confidence_scores <- confidence_scores / max(confidence_scores, na.rm = TRUE)
    }
  } else {
    # Fallback confidence score
    confidence_scores <- rep(0.5, length(cell_predictions))
  }

  # Handle NA values
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0

  # Ensure consistent lengths
  n_cells <- ncol(seurat_test)
  if(length(cell_predictions) != n_cells) {
    cell_predictions <- rep("Unknown", n_cells)
    confidence_scores <- rep(0, n_cells)
  }

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    peak_memory_mb = peak_memory_mb
  ))

  }, error = function(e) {
    warning(paste("singleCellNet error:", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_singleCellNet <- run_singleCellNet_function