# benchmarking_helpers.R
#################################################
# Helper Functions for Cell Type Annotation Benchmarking
# Contains utility functions for metrics, cross-validation, and result aggregation
#################################################

#' Calculate Comprehensive Performance Metrics
#'
#' Purpose: Computes accuracy, precision, recall, F1, macro F1, Cohen's Kappa,
#'          MCC, rare type F1, and unassigned rate
#' Inputs:
#'   - predicted: Vector of predicted cell types
#'   - true_labels: Vector of ground truth cell types
#' Outputs: List with overall_accuracy, macro_f1, cohens_kappa, per_class_metrics,
#'          mcc, rare_type_f1, rare_type_name, unassigned_rate
#' Key Logic: Handles per-class TP/FP/FN calculations, excludes "Unknown" predictions
calculate_metrics <- function(predicted, true_labels) {
  # Overall Accuracy
  overall_accuracy <- mean(predicted == true_labels, na.rm = TRUE)

  # Unassigned rate — computed BEFORE filtering out "Unknown"
  unassigned_rate <- sum(predicted == "Unknown", na.rm = TRUE) / length(predicted)

  # Per-class metrics
  cell_types <- unique(c(predicted, true_labels))
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "Unknown"]
  
  precision_scores <- numeric(length(cell_types))
  recall_scores <- numeric(length(cell_types))
  f1_scores <- numeric(length(cell_types))
  
  for(i in 1:length(cell_types)) {
    ct <- cell_types[i]
    
    # TP, FP, FN calculations
    TP <- sum(predicted == ct & true_labels == ct, na.rm = TRUE)
    FP <- sum(predicted == ct & true_labels != ct, na.rm = TRUE)  
    FN <- sum(predicted != ct & true_labels == ct, na.rm = TRUE)
    
    # Precision and Recall
    precision_scores[i] <- if(TP + FP > 0) TP / (TP + FP) else 0
    recall_scores[i] <- if(TP + FN > 0) TP / (TP + FN) else 0
    
    # F1 Score
    if(precision_scores[i] + recall_scores[i] > 0) {
      f1_scores[i] <- 2 * (precision_scores[i] * recall_scores[i]) / (precision_scores[i] + recall_scores[i])
    } else {
      f1_scores[i] <- 0
    }
  }
  
  # Macro F1 Score (unweighted average)
  macro_f1 <- mean(f1_scores, na.rm = TRUE)

  # Cohen's Kappa - measures agreement adjusted for chance
  # Calculate expected agreement (p_e)
  n <- length(predicted)
  p_e <- 0
  for(ct in cell_types) {
    p_pred <- sum(predicted == ct, na.rm = TRUE) / n
    p_true <- sum(true_labels == ct, na.rm = TRUE) / n
    p_e <- p_e + (p_pred * p_true)
  }

  # Calculate kappa
  cohens_kappa <- if(1 - p_e > 0) {
    (overall_accuracy - p_e) / (1 - p_e)
  } else {
    1  # Perfect agreement if expected agreement is 1
  }

  # MCC - Multiclass Matthews Correlation Coefficient (Gorodkin 2004)
  # Uses confusion matrix C directly; formula:
  #   MCC = (c*N - sum_k t_k*p_k) / sqrt((N^2 - sum_k t_k^2) * (N^2 - sum_k p_k^2))
  valid_mask <- predicted != "Unknown" & !is.na(predicted)
  pred_valid <- predicted[valid_mask]
  true_valid <- true_labels[valid_mask]
  if (length(pred_valid) > 0) {
    all_levels <- union(unique(pred_valid), unique(true_valid))
    C <- as.matrix(table(factor(pred_valid, levels = all_levels),
                         factor(true_valid, levels = all_levels)))
    K <- nrow(C)
    N_mcc <- sum(C)
    c_trace <- sum(diag(C))
    t_k <- rowSums(C)  # predicted counts per class
    p_k <- colSums(C)  # true counts per class
    mcc_num <- c_trace * N_mcc - sum(t_k * p_k)
    mcc_den <- sqrt((N_mcc^2 - sum(t_k^2)) * (N_mcc^2 - sum(p_k^2)))
    mcc <- if (mcc_den > 0) mcc_num / mcc_den else 0
  } else {
    mcc <- 0
  }

  # Rare Type F1 — F1 of the single rarest cell type in true_labels
  per_class_df <- data.frame(
    cell_type = cell_types,
    precision = precision_scores,
    recall = recall_scores,
    f1 = f1_scores
  )
  true_counts <- table(true_labels[true_labels != "Unknown" & !is.na(true_labels)])
  if (length(true_counts) > 0) {
    rare_type_name <- names(which.min(true_counts))
    rare_idx <- which(per_class_df$cell_type == rare_type_name)
    rare_type_f1 <- if (length(rare_idx) > 0) per_class_df$f1[rare_idx[1]] else 0
  } else {
    rare_type_name <- NA
    rare_type_f1 <- 0
  }

  return(list(
    overall_accuracy = overall_accuracy,
    macro_f1 = macro_f1,
    cohens_kappa = cohens_kappa,
    per_class_metrics = per_class_df,
    mcc = mcc,
    rare_type_f1 = rare_type_f1,
    rare_type_name = rare_type_name,
    unassigned_rate = unassigned_rate
  ))
}

#' Create K-Fold Cross-Validation Splits
#' 
#' Purpose: Creates k folds either grouped by donor or stratified by cell type
#' Inputs:
#'   - seurat_obj: Seurat object with metadata
#'   - k: Number of folds (default 5)
#'   - group_by: Column for grouped CV (e.g. "donor_id"), NA for stratified
#'   - stratify_by: Column for stratified CV (e.g. "Ground_Truth_Celltype")
#'   - seed: Random seed for reproducibility
#' Outputs: List of k vectors, each containing cell IDs for test set of that fold
#' Key Logic: Grouped mode keeps donor samples together, stratified mode balances cell types
create_cv_folds <- function(seurat_obj, k = 5, group_by = "donor_id", stratify_by = "Ground_Truth_Celltype", seed = 123) {
  set.seed(seed)
  
  fold_cell_lists <- list()
  
  if (!is.na(group_by) && group_by %in% colnames(seurat_obj@meta.data)) {
    # --- GROUPED MODE ---
    cat(sprintf("Mode: Grouped K-Fold CV by '%s'.\n", group_by))
    
    all_groups <- unique(seurat_obj@meta.data[[group_by]])
    if (length(all_groups) < k) stop("Number of groups is less than k.")
    
    shuffled_groups <- sample(all_groups)
    group_folds <- cut(seq_along(shuffled_groups), breaks = k, labels = FALSE)
    
    for (i in 1:k) {
      test_groups <- shuffled_groups[group_folds == i]
      # Get cell barcodes belonging to the test groups for this fold
      fold_cell_lists[[i]] <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_by]] %in% test_groups, ])
    }
    
  } else {
    # --- STRATIFIED MODE ---
    cat(sprintf("Mode: Stratified K-Fold CV by '%s'.\n", stratify_by))
    
    if (!stratify_by %in% colnames(seurat_obj@meta.data)) {
      stop("Stratification variable not found in metadata.")
    }
    
    # createFolds returns training indices, we need test indices
    all_indices <- 1:ncol(seurat_obj)
    train_indices_list <- createFolds(factor(seurat_obj@meta.data[[stratify_by]]), 
                                     k = k, list = TRUE, returnTrain = TRUE)
    
    for(i in 1:k) {
      test_indices <- setdiff(all_indices, train_indices_list[[i]])
      fold_cell_lists[[i]] <- colnames(seurat_obj)[test_indices]
    }
  }
  
  cat("Fold creation complete.\n")
  
  # Print fold statistics
  for(i in 1:k) {
    cat(sprintf("Fold %d: %d cells\n", i, length(fold_cell_lists[[i]])))
  }
  
  return(fold_cell_lists)
}

# Cache system for fold data
fold_cache <- list()

#' Get Fold Data with Caching
#' 
#' Purpose: Splits data into train/test for a fold and finds marker genes with caching
#' Inputs:
#'   - seurat_obj: Full Seurat object
#'   - folds: Output from create_cv_folds()
#'   - fold_index: Which fold to process (1 to k)
#' Outputs: List with train (Seurat), test (Seurat), markers (dataframe), fold_info
#' Key Logic: Expensive FindAllMarkers() is cached per fold to avoid recomputation
get_fold_data_cached <- function(seurat_obj, folds, fold_index) {
  cache_key <- paste0("fold_", fold_index)
  
  if(!cache_key %in% names(fold_cache)) {
    cat(sprintf("    Computing fold %d data and markers...\n", fold_index))
    
    # Split data
    test_cells <- folds[[fold_index]]
    train_cells <- setdiff(colnames(seurat_obj), test_cells)
    
    seurat_train <- subset(seurat_obj, cells = train_cells)
    seurat_test <- subset(seurat_obj, cells = test_cells)
    
    cat(sprintf("    Training: %d cells, Testing: %d cells\n", 
                ncol(seurat_train), ncol(seurat_test)))
    
    # Find markers (expensive operation - do once per fold)
    Idents(seurat_train) <- seurat_train$Ground_Truth_Celltype
   
     markers <- tryCatch({
      FindAllMarkers(seurat_train,
                     only.pos = FALSE,  # Include both pos/neg markers
                     verbose = FALSE,
                     group.by = "ident",
                     min.cells.group = 3)
    }, error = function(e) {
      warning("FindAllMarkers failed, using empty marker list: ", e$message)
      data.frame()  # Return empty data frame if it fails
    })
    
    cat(sprintf("Found %d marker genes across %d cell types\n", 
                nrow(markers), length(unique(markers$cluster))))
    
    # Cache results
    fold_cache[[cache_key]] <<- list(
      train = seurat_train,
      test = seurat_test, 
      markers = markers,
      fold_info = list(
        train_cells = length(train_cells),
        test_cells = length(test_cells),
        train_cell_types = table(seurat_train$Ground_Truth_Celltype),
        test_cell_types = table(seurat_test$Ground_Truth_Celltype)
      )
    )
  } else {
    cat(sprintf("    Using cached fold %d data\n", fold_index))
  }
  
  return(fold_cache[[cache_key]])
}

#' Aggregate Results from Single Tool Across All Folds
#'
#' Purpose: Combines CV results for one tool into summary statistics
#' Inputs:
#'   - all_predictions: Combined predictions from all folds
#'   - all_true_labels: Combined true labels from all folds
#'   - all_confidence_scores: Combined confidence scores
#'   - all_cell_ids: Combined cell identifiers
#'   - fold_metrics: List of per-fold metric results
#'   - fold_runtimes: Vector of per-fold runtimes
#'   - fold_peak_system_memories: Vector of per-fold peak system memory usage (MB, peak RSS)
#'   - tool_name: Name of the tool
#' Outputs: List with pooled_metrics, fold_variation, runtime_stats, system_memory_stats, detailed_results
#' Key Logic: Provides both pooled analysis (Method B) and fold variation stats (Method A)
aggregate_tool_results <- function(all_predictions, all_true_labels, all_confidence_scores,
                                  all_cell_ids, fold_metrics, fold_runtimes, fold_peak_system_memories, tool_name) {
  
  # Remove failed folds
  valid_folds <- !sapply(fold_metrics, is.null)
  fold_metrics <- fold_metrics[valid_folds]
  fold_runtimes <- fold_runtimes[!is.na(fold_runtimes)]
  
  # Pooled metrics (Method B - preferred)
  pooled_metrics <- calculate_metrics(all_predictions, all_true_labels)
  pooled_metrics$confusion_matrix <- table(True = all_true_labels, Predicted = all_predictions)
  
  # Fold-wise variation (Method A - for stability assessment)
  if(length(fold_metrics) > 0) {
    fold_accuracies <- sapply(fold_metrics, function(x) x$overall_accuracy)
    fold_f1s <- sapply(fold_metrics, function(x) x$macro_f1)
    fold_kappas <- sapply(fold_metrics, function(x) x$cohens_kappa)
    fold_mccs <- sapply(fold_metrics, function(x) x$mcc)
    fold_rare_f1s <- sapply(fold_metrics, function(x) x$rare_type_f1)
    fold_unassigned <- sapply(fold_metrics, function(x) x$unassigned_rate)

    fold_variation <- list(
      accuracy_mean = mean(fold_accuracies, na.rm = TRUE),
      accuracy_std = sd(fold_accuracies, na.rm = TRUE),
      accuracy_folds = fold_accuracies,
      accuracy_ci = mean(fold_accuracies, na.rm = TRUE) + c(-1.96, 1.96) * sd(fold_accuracies, na.rm = TRUE) / sqrt(length(fold_accuracies)),

      f1_mean = mean(fold_f1s, na.rm = TRUE),
      f1_std = sd(fold_f1s, na.rm = TRUE),
      f1_folds = fold_f1s,
      f1_ci = mean(fold_f1s, na.rm = TRUE) + c(-1.96, 1.96) * sd(fold_f1s, na.rm = TRUE) / sqrt(length(fold_f1s)),

      kappa_mean = mean(fold_kappas, na.rm = TRUE),
      kappa_std = sd(fold_kappas, na.rm = TRUE),
      kappa_folds = fold_kappas,
      kappa_ci = mean(fold_kappas, na.rm = TRUE) + c(-1.96, 1.96) * sd(fold_kappas, na.rm = TRUE) / sqrt(length(fold_kappas)),

      mcc_mean = mean(fold_mccs, na.rm = TRUE),
      mcc_std = sd(fold_mccs, na.rm = TRUE),
      mcc_folds = fold_mccs,
      mcc_ci = mean(fold_mccs, na.rm = TRUE) + c(-1.96, 1.96) * sd(fold_mccs, na.rm = TRUE) / sqrt(length(fold_mccs)),

      rare_f1_mean = mean(fold_rare_f1s, na.rm = TRUE),
      rare_f1_std = sd(fold_rare_f1s, na.rm = TRUE),
      rare_f1_folds = fold_rare_f1s,
      rare_f1_ci = mean(fold_rare_f1s, na.rm = TRUE) + c(-1.96, 1.96) * sd(fold_rare_f1s, na.rm = TRUE) / sqrt(length(fold_rare_f1s)),

      unassigned_mean = mean(fold_unassigned, na.rm = TRUE),
      unassigned_std = sd(fold_unassigned, na.rm = TRUE)
    )
  } else {
    fold_variation <- list(
      accuracy_mean = NA, accuracy_std = NA, accuracy_folds = c(),
      f1_mean = NA, f1_std = NA, f1_folds = c(),
      kappa_mean = NA, kappa_std = NA, kappa_folds = c(),
      mcc_mean = NA, mcc_std = NA, mcc_folds = c(), mcc_ci = c(NA, NA),
      rare_f1_mean = NA, rare_f1_std = NA, rare_f1_folds = c(), rare_f1_ci = c(NA, NA),
      unassigned_mean = NA, unassigned_std = NA
    )
  }
  
  # Runtime statistics
  runtime_stats <- list(
    mean_seconds = mean(fold_runtimes, na.rm = TRUE),
    std_seconds = sd(fold_runtimes, na.rm = TRUE),
    total_seconds = sum(fold_runtimes, na.rm = TRUE),
    fold_runtimes = fold_runtimes
  )

  # System memory statistics (peak RSS per fold)
  valid_sys_memories <- fold_peak_system_memories[!is.na(fold_peak_system_memories)]
  system_memory_stats <- list(
    mean_mb            = if (length(valid_sys_memories) > 0) mean(valid_sys_memories) else NA,
    std_mb             = if (length(valid_sys_memories) > 1) sd(valid_sys_memories) else NA,
    max_mb             = if (length(valid_sys_memories) > 0) max(valid_sys_memories) else NA,
    fold_peak_memories = valid_sys_memories
  )

  return(list(
    tool_name = tool_name,
    pooled_metrics = pooled_metrics,
    fold_variation = fold_variation,
    runtime_stats = runtime_stats,
    system_memory_stats = system_memory_stats,
    detailed_results = list(
      predictions = all_predictions,
      true_labels = all_true_labels,
      confidence_scores = all_confidence_scores,
      cell_ids = all_cell_ids
    ),
    num_successful_folds = length(fold_metrics),
    num_total_folds = length(valid_folds)
  ))
}

#' Run Cross-Validation for Single Tool
#' 
#' Purpose: Executes one tool across all CV folds with error handling
#' Inputs:
#'   - seurat_obj: Full Seurat object
#'   - tool_function: Function that takes (train, test, markers) and returns predictions
#'   - folds: Output from create_cv_folds()
#'   - tool_name: Name for logging/identification
#' Outputs: Aggregated results from aggregate_tool_results()
#' Key Logic: Runs tool on each fold, accumulates results, handles failures gracefully
run_single_tool_cv <- function(seurat_obj, tool_function, folds, tool_name) {
  cat(sprintf("\n=== Running %s ===\n", tool_name))
  
  # Storage for this tool
  all_predictions <- c()
  all_true_labels <- c()
  all_confidence_scores <- c()
  all_cell_ids <- c()
  fold_metrics <- list()
  fold_runtimes <- c()
  fold_peak_system_memories <- c()
  
  # Process each fold
  for(i in 1:length(folds)) {
    cat(sprintf("  Fold %d/%d...\n", i, length(folds)))
    
    # Get pre-computed fold data
    fold_data <- get_fold_data_cached(seurat_obj, folds, i)
    
    # Run tool on this fold
    start_time <- Sys.time()
    tryCatch({
      fold_result <- tool_function(fold_data$train, fold_data$test, fold_data$markers)
      runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Validate tool output
      if(!all(c("predictions", "true_labels", "confidence_scores", "cell_ids") %in% names(fold_result))) {
        stop("Tool function must return list with: predictions, true_labels, confidence_scores, cell_ids")
      }
      
      # Accumulate results for pooled analysis
      all_predictions <- c(all_predictions, fold_result$predictions)
      all_true_labels <- c(all_true_labels, fold_result$true_labels)
      all_confidence_scores <- c(all_confidence_scores, fold_result$confidence_scores)
      all_cell_ids <- c(all_cell_ids, fold_result$cell_ids)
      
      # Store fold-wise metrics
      fold_metrics[[i]] <- calculate_metrics(fold_result$predictions, fold_result$true_labels)

      # Use tool-reported runtime if available, fall back to outer wall clock
      tool_runtime <- fold_result$runtime_secs
      if (is.null(tool_runtime) || is.na(tool_runtime)) tool_runtime <- runtime
      fold_runtimes[i] <- tool_runtime

      # Extract peak system memory if available
      peak_sys_mem <- fold_result$peak_system_memory_mb
      if (is.null(peak_sys_mem)) peak_sys_mem <- NA
      fold_peak_system_memories[i] <- peak_sys_mem

      cat(sprintf("    Completed in %.2f seconds, accuracy: %.2f%%, memory: %.1f MB\n",
                  tool_runtime, fold_metrics[[i]]$overall_accuracy * 100,
                  ifelse(is.na(peak_sys_mem), 0, peak_sys_mem)))
      
    }, error = function(e) {
      warning(sprintf("Tool %s failed on fold %d: %s", tool_name, i, e$message))
      fold_metrics[[i]] <- NULL
      fold_runtimes[i] <- NA
      fold_peak_system_memories[i] <- NA
      cat(sprintf("    FAILED: %s\n", e$message))
    })
  }

  # Aggregate results for this tool
  return(aggregate_tool_results(all_predictions, all_true_labels,
                               all_confidence_scores, all_cell_ids,
                               fold_metrics, fold_runtimes, fold_peak_system_memories, tool_name))
}

#' Aggregate Results from All Tools 
#' 
#' Purpose: Creates comparison table and summary statistics across all tools
#' Inputs:
#'   - all_tool_results: List of results from run_single_tool_cv() for each tool
#' Outputs: List with comparison_table, detailed_results, best_tool, summary_stats
#' Key Logic: Extracts key metrics from each tool, sorts by accuracy, provides rankings
aggregate_all_tools <- function(all_tool_results) {
  
  # Create comparison table
  comparison_df <- data.frame(
    tool = character(),
    pooled_accuracy = numeric(),
    pooled_macro_f1 = numeric(),
    pooled_kappa = numeric(),
    pooled_mcc = numeric(),
    pooled_rare_type_f1 = numeric(),
    pooled_unassigned_rate = numeric(),
    accuracy_mean = numeric(),
    accuracy_std = numeric(),
    accuracy_ci_lower = numeric(),
    accuracy_ci_upper = numeric(),
    f1_mean = numeric(),
    f1_std = numeric(),
    f1_ci_lower = numeric(),
    f1_ci_upper = numeric(),
    kappa_mean = numeric(),
    kappa_std = numeric(),
    kappa_ci_lower = numeric(),
    kappa_ci_upper = numeric(),
    mcc_mean = numeric(),
    mcc_std = numeric(),
    rare_f1_mean = numeric(),
    rare_f1_std = numeric(),
    unassigned_mean = numeric(),
    unassigned_std = numeric(),
    runtime_mean = numeric(),
    runtime_std = numeric(),
    peak_system_memory_mean_mb = numeric(),
    peak_system_memory_std_mb = numeric(),
    peak_system_memory_max_mb = numeric(),
    successful_folds = integer(),
    stringsAsFactors = FALSE
  )
  
  # Extract metrics for each tool
  for(tool_name in names(all_tool_results)) {
    result <- all_tool_results[[tool_name]]
    
    comparison_df <- rbind(comparison_df, data.frame(
      tool = tool_name,
      pooled_accuracy = result$pooled_metrics$overall_accuracy,
      pooled_macro_f1 = result$pooled_metrics$macro_f1,
      pooled_kappa = result$pooled_metrics$cohens_kappa,
      pooled_mcc = result$pooled_metrics$mcc,
      pooled_rare_type_f1 = result$pooled_metrics$rare_type_f1,
      pooled_unassigned_rate = result$pooled_metrics$unassigned_rate,
      accuracy_mean = result$fold_variation$accuracy_mean,
      accuracy_std = result$fold_variation$accuracy_std,
      accuracy_ci_lower = result$fold_variation$accuracy_ci[1],
      accuracy_ci_upper = result$fold_variation$accuracy_ci[2],
      f1_mean = result$fold_variation$f1_mean,
      f1_std = result$fold_variation$f1_std,
      f1_ci_lower = result$fold_variation$f1_ci[1],
      f1_ci_upper = result$fold_variation$f1_ci[2],
      kappa_mean = result$fold_variation$kappa_mean,
      kappa_std = result$fold_variation$kappa_std,
      kappa_ci_lower = result$fold_variation$kappa_ci[1],
      kappa_ci_upper = result$fold_variation$kappa_ci[2],
      mcc_mean = result$fold_variation$mcc_mean,
      mcc_std = result$fold_variation$mcc_std,
      rare_f1_mean = result$fold_variation$rare_f1_mean,
      rare_f1_std = result$fold_variation$rare_f1_std,
      unassigned_mean = result$fold_variation$unassigned_mean,
      unassigned_std = result$fold_variation$unassigned_std,
      runtime_mean = result$runtime_stats$mean_seconds,
      runtime_std = result$runtime_stats$std_seconds,
      peak_system_memory_mean_mb = result$system_memory_stats$mean_mb,
      peak_system_memory_std_mb = result$system_memory_stats$std_mb,
      peak_system_memory_max_mb = result$system_memory_stats$max_mb,
      successful_folds = result$num_successful_folds
    ))
  }
  
  # Sort by pooled kappa (primary metric)
  comparison_df <- comparison_df[order(-comparison_df$pooled_kappa), ]

  return(list(
    comparison_table = comparison_df,
    detailed_results = all_tool_results,
    best_tool = comparison_df$tool[1],
    summary_stats = list(
      num_tools = nrow(comparison_df),
      best_accuracy = max(comparison_df$pooled_accuracy, na.rm = TRUE),
      best_f1 = max(comparison_df$pooled_macro_f1, na.rm = TRUE),
      best_kappa = max(comparison_df$pooled_kappa, na.rm = TRUE),
      best_mcc = max(comparison_df$pooled_mcc, na.rm = TRUE),
      fastest_tool = comparison_df$tool[which.min(comparison_df$runtime_mean)]
    )
  ))
}

cat("Benchmarking helper functions loaded successfully!\n")