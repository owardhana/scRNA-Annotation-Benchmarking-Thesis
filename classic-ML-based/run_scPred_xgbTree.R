# run_scPred_xgbTree.R
#################################################
# scPred Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' scPred Cell Type Annotation Function (xgbTree)
#'
#' Purpose: Run scPred algorithm using xgboost-based classification
#' Inputs:
#'   - seurat_train: Training Seurat object (used to train ML model)
#'   - seurat_test: Test Seurat object to predict
#'   - markers: Marker genes dataframe from FindAllMarkers() (for interface consistency)
#' Outputs: List with predictions, true_labels, confidence_scores, cell_ids
#' Algorithm: Uses scPred's getFeatureSpace → trainModel → scPredict workflow.
#'            Uses a custom caret model definition to bypass the broken caret xgbTree
#'            interface (ALTREP/DMatrix incompatibility). Memory optimised via reduced
#'            grid search, immediate DMatrix cleanup, and nthread = 1.
run_scPred_xgbTree_function <- function(seurat_train, seurat_test, markers) {

  # Load required libraries
  if (!requireNamespace("scPred", quietly = TRUE)) {
    stop("scPred package not available. Please install it first.")
  }
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost package not available. Please install it first.")
  }
  library(scPred)
  library(Seurat)
  library(caret)
  library(xgboost)

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

  # Track peak memory usage
  runtime_secs <- NA
  peak_system_memory_mb <- NA

  # Custom caret model definition — bypasses broken caret xgbTree interface.
  # Grid is kept small (4 combos) to minimise caret's internal CV memory footprint.
  customXgb <- list(
    library = "xgboost",
    type = "Classification",
    parameters = data.frame(
      parameter = c("nrounds", "max_depth", "eta", "gamma",
                    "colsample_bytree", "min_child_weight", "subsample"),
      class = rep("numeric", 7),
      label = c("# Boosting Iterations", "Max Tree Depth", "Shrinkage",
                "Minimum Loss Reduction", "Subsample Ratio of Columns",
                "Minimum Sum of Instance Weight", "Subsample Percentage")
    ),
    grid = function(x, y, len = NULL, search = "grid") {
      expand.grid(
        nrounds          = c(100, 200),
        max_depth        = c(3, 6),
        eta              = 0.3,
        gamma            = 0,
        colsample_bytree = 0.8,
        min_child_weight = 1,
        subsample        = 0.8
      )
    },
    loop = NULL,
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      x <- as.matrix(x)
      storage.mode(x) <- "double"

      # Ensure valid column names
      if (is.null(colnames(x))) {
        colnames(x) <- paste0("V", seq_len(ncol(x)))
      }
      feat_names <- make.names(colnames(x))
      colnames(x) <- feat_names

      # Prepare labels
      if (is.factor(y)) {
        y_numeric <- as.integer(y) - 1L
        num_class  <- length(levels(y))
        lev_names  <- levels(y)
      } else {
        y_numeric <- y
        num_class  <- NULL
        lev_names  <- NULL
      }

      # Build xgboost parameter list
      # nthread = 1: prevents per-fold thread explosion when called inside caret CV
      params <- list(
        eta              = param$eta,
        max_depth        = param$max_depth,
        gamma            = param$gamma,
        colsample_bytree = param$colsample_bytree,
        min_child_weight = param$min_child_weight,
        subsample        = param$subsample,
        tree_method      = "hist",
        nthread          = 1
      )

      if (!is.null(num_class) && num_class > 2L) {
        params$objective   <- "multi:softprob"
        params$num_class   <- num_class
        params$eval_metric <- "mlogloss"
      } else if (!is.null(num_class)) {
        params$objective   <- "binary:logistic"
        params$eval_metric <- "logloss"
      } else {
        params$objective   <- "reg:squarederror"
        params$eval_metric <- "rmse"
      }

      # Create DMatrix then immediately free the dense matrix
      dtrain <- xgboost::xgb.DMatrix(data = x, label = y_numeric)
      rm(x); gc()

      model <- xgboost::xgb.train(
        params  = params,
        data    = dtrain,
        nrounds = param$nrounds,
        verbose = 0
      )
      rm(dtrain); gc()

      model_wrapper <- list(
        model_obj  = model,
        feat_names = feat_names,
        n_features = length(feat_names),
        levels     = lev_names,
        num_class  = num_class
      )
      class(model_wrapper) <- c("xgb_wrapper", "list")
      return(model_wrapper)
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      newdata <- as.matrix(newdata)
      storage.mode(newdata) <- "double"
      if (ncol(newdata) == modelFit$n_features) {
        colnames(newdata) <- modelFit$feat_names
      }
      dtest <- xgboost::xgb.DMatrix(data = newdata)
      pred_probs <- predict(modelFit$model_obj, dtest)
      rm(dtest)

      if (modelFit$num_class == 2L) {
        pred <- ifelse(pred_probs > 0.5, modelFit$levels[2], modelFit$levels[1])
      } else {
        pred_matrix <- matrix(pred_probs, ncol = modelFit$num_class, byrow = TRUE)
        pred <- modelFit$levels[max.col(pred_matrix)]
      }
      return(pred)
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      newdata <- as.matrix(newdata)
      storage.mode(newdata) <- "double"
      if (ncol(newdata) == modelFit$n_features) {
        colnames(newdata) <- modelFit$feat_names
      }
      dtest <- xgboost::xgb.DMatrix(data = newdata)
      pred_probs <- predict(modelFit$model_obj, dtest)
      rm(dtest)

      if (modelFit$num_class == 2L) {
        prob_matrix <- cbind(1 - pred_probs, pred_probs)
      } else {
        prob_matrix <- matrix(pred_probs, ncol = modelFit$num_class, byrow = TRUE)
      }
      colnames(prob_matrix) <- modelFit$levels
      return(as.data.frame(prob_matrix))
    },
    predictors = function(x, ...) x$feat_names,
    tags = c("Tree-Based Model", "Ensemble Model", "Boosting", "Implicit Feature Selection"),
    levels = function(x) x$levels,
    sort = function(x) x[order(x$nrounds, x$max_depth, x$eta), ]
  )

  # Run scPred workflow with memory tracking
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    warning("peakRAM package not available for memory/time tracking")

      start_time <- Sys.time()
    # 1. Create feature space from training data
    seurat_train <- getFeatureSpace(seurat_train, "Ground_Truth_Celltype")

    # 2. Train model — allowParallel = FALSE prevents per-worker data copies
    seurat_train <- trainModel(seurat_train, model = customXgb, allowParallel = FALSE)

    # 3. Predict on test data
    seurat_test <- scPredict(seurat_test, seurat_train)
    runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  } else {
    library(peakRAM)
    peakRAM_result <- peakRAM::peakRAM(
      {
        # 1. Create feature space from training data
        seurat_train <- getFeatureSpace(seurat_train, "Ground_Truth_Celltype")

        # 2. Train model — allowParallel = FALSE prevents per-worker data copies
        seurat_train <- trainModel(seurat_train, model = customXgb, allowParallel = FALSE)

        # 3. Predict on test data
        seurat_test <- scPredict(seurat_test, seurat_train)
      
    })
    runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
    peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
  }

  # Extract predictions and confidence scores
  cell_predictions <- seurat_test$scpred_prediction
  true_labels <- seurat_test$Ground_Truth_Celltype

  # Extract confidence scores from probabilities stored in metadata
  prob_cols <- grep("scpred_", colnames(seurat_test@meta.data), value = TRUE)
  prob_cols <- prob_cols[prob_cols != "scpred_prediction"]

  if (length(prob_cols) > 0) {
    prob_matrix <- seurat_test@meta.data[, prob_cols, drop = FALSE]
    confidence_scores <- apply(prob_matrix, 1, max, na.rm = TRUE)
  } else {
    confidence_scores <- rep(1.0, length(cell_predictions))
    print("No confidence found")
    print(seurat_test@meta.data)
  }

  # Handle NA values
  cell_predictions[is.na(cell_predictions)] <- "Unknown"
  confidence_scores[is.na(confidence_scores)] <- 0

  # Return standardized format
  return(list(
    predictions = as.character(cell_predictions),
    true_labels = true_labels,
    confidence_scores = confidence_scores,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))

  }, error = function(e) {
    warning(paste("scPred (xgbTree) error:", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_scPred_xgbTree <- run_scPred_xgbTree_function
