# Main_Benchmarking.R
#################################################
# Scalable Benchmarking Framework for Cell Type Annotation Tools
# 
# This is the main execution file that orchestrates the benchmarking process:
# 1. Loads configuration and tool registry
# 2. Sources helper functions from benchmarking_helpers.R
# 3. Provides main execution function run_all_tools_cv()
# 4. Contains example usage patterns
#
# Supports both marker-based and reference-based methods
# Uses k-fold cross-validation with comprehensive result aggregation
#################################################
# Load required libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(caret)
library(SingleCellExperiment)
library(SummarizedExperiment)
#library(monocle)
library(tensorflow)
library(reticulate)
library(bench)
library(pryr)



# Load helper functions
source("benchmarking_helpers.R")

# Load tool functions
source("marker-based/run_SCINA.R")
source("marker-based/run_scSorter.R")
source("marker-based/run_scType.R")
source("marker-based/run_scCATCH.R")
source("marker-based/run_Garnett.R")
source("marker-based/run_SCSA.R")
source("marker-based/run_clustifyr_hyper.R")
source("marker-based/run_clustifyr_jaccard.R")
source("marker-based/run_ScInfeR.R")

source("reference-based/run_SingleR.R")
source("reference-based/run_scibetR.R")
source("reference-based/run_scmap_cell.R")
source("reference-based/run_scmap_cluster.R")
source("reference-based/run_Seurat_Transfer_PCA.R")
source("reference-based/run_Seurat_Transfer_RPCA.R")
source("reference-based/run_Seurat_Transfer_CCA.R")
source("reference-based/run_CIPR.R")

source("classic-ML-based/run_scPred.R")
source("classic-ML-based/run_scPred_avNNet.R")
source("classic-ML-based/run_scPred_xgbTree.R")
source("classic-ML-based/run_scPred_rf.R")
source("classic-ML-based/run_scPred_glm.R")
source("classic-ML-based/run_scPred_glmboost.R")
source("classic-ML-based/run_scPred_adaboost.R")
source("classic-ML-based/run_scPred_lda.R")
source("classic-ML-based/run_scPred_knn.R")
source("classic-ML-based/run_CHETAH.R")
source("classic-ML-based/run_scClassify.R")
source("classic-ML-based/run_scAnnotatR.R") 
source("classic-ML-based/run_singleCellNet.R")
source("classic-ML-based/run_CellTypist.R")
source("classic-ML-based/run_scAnnotate.R")
source("classic-ML-based/run_CALLR.R")
source("classic-ML-based/run_scAnno.R")
source("classic-ML-based/run_CaSTLe.R")
source("classic-ML-based/run_scID.R")
source("classic-ML-Based/run_scPred_svmLinear.R")
source("classic-ML-Based/run_scPred_nb.R")
source("classic-ML-Based/run_scPred_svmPoly.R")
source("classic-ML-Based/run_scPred_bayesglm.R")
source("classic-ML-Based/run_scPred_earth.R")
source("classic-ML-Based/run_scPred_mlp.R")
source("classic-ML-Based/run_scPred_nnet.R")
source("classic-ML-Based/run_scPred_regLogistic.R")
source("classic-ML-Based/run_scPred_multinom.R")
source("classic-ML-Based/run_scPred_glmnet.R")

source("DL-based/run_scLearn.R")
source("DL-based/run_CAMLU.R")
source("DL-based/run_NeuCA.R")
source("DL-based/run_scPred_mxnetAdam.R")

source("LLM-based/run_mLLMCelltype.R")
source("LLM-based/run_CASSIA.R")
source("LLM-based/run_GPTCelltype.R")

source("LLM-based/run_mLLMCelltype_claude_sonnet4.5.R")
source("LLM-based/run_mLLMCelltype_gpt5.R")
source("LLM-based/run_mLLMCelltype_gemini_pro2.5.R")
source("LLM-based/run_mLLMCelltype_grok4_fast.R")
source("LLM-based/run_mLLMCelltype_deepseek_v3.1_terminus.R")
source("LLM-based/run_mLLMCelltype_qwen3_max.R")
source("LLM-based/run_mLLMCelltype_llama4_maverick.R")


#################################################
# TOOL REGISTRY AND EXECUTION
#################################################

# Tool registry - will be populated as we convert each tool to function format
TOOL_REGISTRY <- list()

# Add converted tool functions
TOOL_REGISTRY[["SCINA"]] <- run_SCINA_function
TOOL_REGISTRY[["scSorter"]] <- run_scSorter_function
TOOL_REGISTRY[["scType"]] <- run_scType_function
TOOL_REGISTRY[["scCATCH"]] <- run_scCATCH_function
TOOL_REGISTRY[["SCSA"]] <- run_SCSA_function
TOOL_REGISTRY[["Garnett"]] <- run_Garnett_function
TOOL_REGISTRY[["clustifyr_hyper"]] <- run_clustifyr_hyper_function
TOOL_REGISTRY[["clustifyr_jaccard"]] <- run_clustifyr_jaccard_function
TOOL_REGISTRY[["ScInfeR"]] <- run_ScInfeR_function

TOOL_REGISTRY[["SingleR"]] <- run_SingleR_function
TOOL_REGISTRY[["Seurat_Transfer_PCA"]] <- run_Seurat_Transfer_PCA_function
TOOL_REGISTRY[["Seurat_Transfer_RPCA"]] <- run_Seurat_Transfer_RPCA_function
TOOL_REGISTRY[["Seurat_Transfer_CCA"]] <- run_Seurat_Transfer_CCA_function
TOOL_REGISTRY[["scmap_cell"]] <- run_scmap_cell_function
TOOL_REGISTRY[["scmap_cluster"]] <- run_scmap_cluster_function
TOOL_REGISTRY[["scibetR"]] <- run_scibetR_function
TOOL_REGISTRY[["CIPR"]] <- run_CIPR_function

TOOL_REGISTRY[["scPred"]] <- run_scPred_function
TOOL_REGISTRY[["scPred_avNNet"]] <- run_scPred_avNNet_function
TOOL_REGISTRY[["scPred_xgbTree"]] <- run_scPred_xgbTree_function
TOOL_REGISTRY[["scPred_rf"]] <- run_scPred_rf_function
TOOL_REGISTRY[["scPred_glm"]] <- run_scPred_glm_function
TOOL_REGISTRY[["scPred_glmboost"]] <- run_scPred_glmboost_function
TOOL_REGISTRY[["scPred_adaboost"]] <- run_scPred_adaboost_function
TOOL_REGISTRY[["scPred_lda"]] <- run_scPred_lda_function
TOOL_REGISTRY[["scPred_knn"]] <- run_scPred_knn_function
TOOL_REGISTRY[["CHETAH"]] <- run_CHETAH_function
TOOL_REGISTRY[["scClassify"]] <- run_scClassify_function
TOOL_REGISTRY[["scAnnotatR"]] <- run_scAnnotatR_function
TOOL_REGISTRY[["singleCellNet"]] <- run_singleCellNet_function
TOOL_REGISTRY[["CellTypist"]] <- run_CellTypist_function
TOOL_REGISTRY[["scAnnotate"]] <- run_scAnnotate_function
TOOL_REGISTRY[["CALLR"]] <- run_CALLR_function
TOOL_REGISTRY[["CaSTLe"]] <- run_CaSTLe_function
TOOL_REGISTRY[["scID"]] <- run_scID_function
TOOL_REGISTRY[["scPred_svmLinear"]]   <- run_scPred_svmLinear_function
TOOL_REGISTRY[["scPred_nb"]]          <- run_scPred_nb_function
TOOL_REGISTRY[["scPred_svmPoly"]]     <- run_scPred_svmPoly_function
TOOL_REGISTRY[["scPred_bayesglm"]]    <- run_scPred_bayesglm_function
TOOL_REGISTRY[["scPred_earth"]]       <- run_scPred_earth_function
TOOL_REGISTRY[["scPred_mlp"]]         <- run_scPred_mlp_function
TOOL_REGISTRY[["scPred_nnet"]]        <- run_scPred_nnet_function
TOOL_REGISTRY[["scPred_regLogistic"]] <- run_scPred_regLogistic_function
TOOL_REGISTRY[["scPred_multinom"]]    <- run_scPred_multinom_function
TOOL_REGISTRY[["scPred_glmnet"]]     <- run_scPred_glmnet_function

TOOL_REGISTRY[["scLearn"]] <- run_scLearn_function
TOOL_REGISTRY[["CAMLU"]] <- run_CAMLU_function
TOOL_REGISTRY[["NeuCA"]] <- run_NeuCA_function
TOOL_REGISTRY[["scPred_mxnetAdam"]]   <- run_scPred_mxnetAdam_function

TOOL_REGISTRY[["mLLMCelltype"]] <- run_mLLMCelltype_function
TOOL_REGISTRY[["CASSIA"]] <- run_CASSIA_function
TOOL_REGISTRY[["GPTCelltype"]] <- run_GPTCelltype_function
TOOL_REGISTRY[["mLLMCelltype_claude_sonnet4.5"]] <- run_mLLMCelltype_claude_sonnet4.5_function
TOOL_REGISTRY[["mLLMCelltype_gpt5"]] <- run_mLLMCelltype_gpt5_function
TOOL_REGISTRY[["mLLMCelltype_gemini_pro2.5"]] <- run_mLLMCelltype_gemini_pro2.5_function
TOOL_REGISTRY[["mLLMCelltype_grok4_fast"]] <- run_mLLMCelltype_grok4_fast_function
TOOL_REGISTRY[["mLLMCelltype_deepseek_v3.1_terminus"]] <- run_mLLMCelltype_deepseek_v3.1_terminus_function
TOOL_REGISTRY[["mLLMCelltype_qwen3_max"]] <- run_mLLMCelltype_qwen3_max_function
TOOL_REGISTRY[["mLLMCelltype_llama4_maverick"]] <- run_mLLMCelltype_llama4_maverick_function



#' Main Benchmarking Execution Function
#' 
#' Purpose: Orchestrates the complete benchmarking workflow
#' Inputs:
#'   - seurat_obj: Seurat object with Ground_Truth_Celltype metadata
#'   - tools_to_run: Vector of tool names to benchmark (default: TOOLS_TO_RUN)
#'   - cv_config: List with k, group_by, stratify_by, seed (default: CV_CONFIG)
#' Outputs: List of results from each tool (pass to aggregate_all_tools())
#' Workflow:
#'   1. Creates k-fold splits once (reused for all tools)
#'   2. Runs each tool across all folds with error handling
#'   3. Returns aggregated results ready for comparison
run_all_tools_cv <- function(seurat_obj, tools_to_run = TOOLS_TO_RUN, 
                             cv_config = CV_CONFIG) {
  
  cat("=== CELL TYPE ANNOTATION BENCHMARKING ===\n")
  cat(sprintf("Dataset: %d cells x %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))
  cat(sprintf("Cell types: %s\n", paste(unique(seurat_obj$Ground_Truth_Celltype), collapse = ", ")))
  cat(sprintf("Tools to run: %s\n", paste(tools_to_run, collapse = ", ")))
  
  # Create folds once, reuse for all tools
  folds <- create_cv_folds(seurat_obj, 
                          k = cv_config$k,
                          group_by = cv_config$group_by, 
                          stratify_by = cv_config$stratify_by,
                          seed = cv_config$seed)
  
  # Initialize results storage
  all_tool_results <- list()
  
  # Run each tool
  for(tool_name in tools_to_run) {
    tool_function <- TOOL_REGISTRY[[tool_name]]
    if(is.null(tool_function)) {
      warning(sprintf("Tool %s not found in registry", tool_name))
      next
    }
    
    # Run CV for this tool
    tool_results <- run_single_tool_cv(seurat_obj, tool_function, folds, tool_name)
    all_tool_results[[tool_name]] <- tool_results
  }
  
  return(all_tool_results)
}

#' Single 80/20 Split Execution Mode
#'
#' Purpose: Runs all tools on a single stratified 80/20 train/test split
#'          instead of k-fold CV. Designed for L9 simulated datasets where
#'          each dataset is an independent Splatter replicate.
#' Inputs:
#'   - seurat_obj: Seurat object with metadata column specified by label_col
#'   - tools_to_run: Vector of tool names to benchmark
#'   - label_col: Metadata column with cell type labels (default: "Ground_Truth_Celltype")
#'   - seed: Random seed for reproducibility
#' Outputs: Named list (tool_name -> result) compatible with aggregate_all_tools()
#'          Each result matches the structure returned by aggregate_tool_results()
run_all_tools_single_split <- function(seurat_obj, tools_to_run = TOOLS_TO_RUN,
                                       label_col = "Ground_Truth_Celltype",
                                       seed = 123) {

  cat("=== CELL TYPE ANNOTATION BENCHMARKING (Single 80/20 Split) ===\n")
  cat(sprintf("Dataset: %d cells x %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))
  cat(sprintf("Cell types: %s\n", paste(unique(seurat_obj@meta.data[[label_col]]), collapse = ", ")))
  cat(sprintf("Tools to run: %s\n", paste(tools_to_run, collapse = ", ")))

  # Create single stratified 80/20 split
  set.seed(seed)
  train_idx <- caret::createDataPartition(seurat_obj@meta.data[[label_col]],
                                          p = 0.80, list = FALSE)
  train_cells <- colnames(seurat_obj)[train_idx]
  test_cells <- colnames(seurat_obj)[-train_idx]

  seurat_train <- subset(seurat_obj, cells = train_cells)
  seurat_test <- subset(seurat_obj, cells = test_cells)

  cat(sprintf("Training: %d cells, Testing: %d cells\n",
              ncol(seurat_train), ncol(seurat_test)))

  # Find markers on training data
  Idents(seurat_train) <- seurat_train@meta.data[[label_col]]
  markers <- tryCatch({
    FindAllMarkers(seurat_train,
                   only.pos = FALSE,
                   verbose = FALSE,
                   group.by = "ident",
                   min.cells.group = 3)
  }, error = function(e) {
    warning("FindAllMarkers failed, using empty marker list: ", e$message)
    data.frame()
  })

  cat(sprintf("Found %d marker genes across %d cell types\n",
              nrow(markers), length(unique(markers$cluster))))

  # Run each tool
  all_tool_results <- list()

  for (tool_name in tools_to_run) {
    tool_function <- TOOL_REGISTRY[[tool_name]]
    if (is.null(tool_function)) {
      warning(sprintf("Tool %s not found in registry", tool_name))
      next
    }

    cat(sprintf("\n=== Running %s ===\n", tool_name))
    start_time <- Sys.time()

    tryCatch({
      tool_result <- tool_function(seurat_train, seurat_test, markers)
      runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      # Validate tool output
      if (!all(c("predictions", "true_labels", "confidence_scores", "cell_ids") %in% names(tool_result))) {
        stop("Tool function must return list with: predictions, true_labels, confidence_scores, cell_ids")
      }

      # Use tool-reported runtime if available
      tool_runtime <- tool_result$runtime_secs
      if (is.null(tool_runtime) || is.na(tool_runtime)) tool_runtime <- runtime

      # Extract peak system memory if available
      peak_sys_mem <- tool_result$peak_system_memory_mb
      if (is.null(peak_sys_mem)) peak_sys_mem <- NA

      # Compute metrics
      metrics <- calculate_metrics(tool_result$predictions, tool_result$true_labels)

      cat(sprintf("    Completed in %.2f seconds, accuracy: %.2f%%, memory: %.1f MB\n",
                  tool_runtime, metrics$overall_accuracy * 100,
                  ifelse(is.na(peak_sys_mem), 0, peak_sys_mem)))

      # Build result structure matching aggregate_tool_results() output
      # Single split: fold_variation has single values, std = NA, ci = c(NA, NA)
      all_tool_results[[tool_name]] <- list(
        tool_name = tool_name,
        pooled_metrics = metrics,
        fold_variation = list(
          accuracy_mean = metrics$overall_accuracy,
          accuracy_std = NA,
          accuracy_folds = metrics$overall_accuracy,
          accuracy_ci = c(NA, NA),
          f1_mean = metrics$macro_f1,
          f1_std = NA,
          f1_folds = metrics$macro_f1,
          f1_ci = c(NA, NA),
          kappa_mean = metrics$cohens_kappa,
          kappa_std = NA,
          kappa_folds = metrics$cohens_kappa,
          kappa_ci = c(NA, NA),
          mcc_mean = metrics$mcc,
          mcc_std = NA,
          mcc_folds = metrics$mcc,
          mcc_ci = c(NA, NA),
          rare_f1_mean = metrics$rare_type_f1,
          rare_f1_std = NA,
          rare_f1_folds = metrics$rare_type_f1,
          rare_f1_ci = c(NA, NA),
          unassigned_mean = metrics$unassigned_rate,
          unassigned_std = NA
        ),
        runtime_stats = list(
          mean_seconds = tool_runtime,
          std_seconds = NA,
          total_seconds = tool_runtime,
          fold_runtimes = tool_runtime
        ),
        system_memory_stats = list(
          mean_mb = peak_sys_mem,
          std_mb = NA,
          max_mb = peak_sys_mem,
          fold_peak_memories = peak_sys_mem
        ),
        detailed_results = list(
          predictions = tool_result$predictions,
          true_labels = tool_result$true_labels,
          confidence_scores = tool_result$confidence_scores,
          cell_ids = tool_result$cell_ids
        ),
        num_successful_folds = 1,
        num_total_folds = 1
      )
      # Add confusion matrix to pooled_metrics
      all_tool_results[[tool_name]]$pooled_metrics$confusion_matrix <-
        table(True = tool_result$true_labels, Predicted = tool_result$predictions)

    }, error = function(e) {
      warning(sprintf("Tool %s failed: %s", tool_name, e$message))
      cat(sprintf("    FAILED: %s\n", e$message))
    })
  }

  return(all_tool_results)
}

#################################################
# CONFIGURATION SECTION
#################################################

# Data paths - Datasets (Seurat with Normalized and Find Variable Features run)
SEURAT_PATH <- "./data/S3_rep2_seurat_subset.rds"

# Cross-validation configuration
CV_CONFIG <- list(
  k = 3,                              # Number of folds
  group_by = NA,              # For grouped CV, set to NA for stratified
  stratify_by = "Ground_Truth_Celltype", # For stratified CV
  seed = 123                          # Reproducibility
)

# L9 simulation split configuration
# Used by run_all_tools_single_split() for simulated datasets
# Each L9 replicate is an independent dataset — no CV needed
SPLIT_CONFIG <- list(
  label_col = "Ground_Truth_Celltype",
  train_prop = 0.80,
  seed = 123
)


# BASIC WORKFLOW:
# 
# 1. Load data
seurat_obj <- readRDS(SEURAT_PATH)
Idents(seurat_obj) <- "Ground_Truth_Celltype"
# 
# 2. Run benchmarking (start with subset for testing)
marker_based <- c(
  "SCINA",
  "scSorter",
  "scType",
  "scCATCH",
  "SCSA",
  #"Garnett", #broken
  "clustifyr_hyper",
  "clustifyr_jaccard",
  "ScInfeR"
)

reference_based <- c(
  "SingleR",
  "Seurat_Transfer_PCA",
  "Seurat_Transfer_RPCA",
  "Seurat_Transfer_CCA",
  "scmap_cell",
  "scmap_cluster",
  "scibetR",
  "CIPR"
)

classicML_based1 <- c(
  "scPred",
  "scPred_avNNet",
  "scPred_xgbTree",
  "scPred_rf",
  "scPred_glm",
  "scPred_glmboost",
  "scPred_lda",
  "scPred_knn",
  "CHETAH",
  "scClassify",
  "scAnnotatR",
  "singleCellNet",
  "CellTypist"
)

classicML_based2 <- c(
  "scAnnotate",
  "CALLR",
  "CaSTLe",
  "scPred_adaboost", 
  "scID"
)

classicML_based3 <- c(
  "scPred_svmLinear",
  "scPred_nb",
  "scPred_svmPoly",
  "scPred_bayesglm",
  "scPred_earth",
  "scPred_mlp",
  "scPred_nnet",
  "scPred_regLogistic",
  "scPred_multinom"
  #"scPred_glmnet" #broken
)

DL_based1 <- c(
  "scLearn",
  "CAMLU",
  "NeuCA"
  #"scPred_mxnetAdam" #broken
)

LLM_Tools <- c(
  "mLLMCelltype",
  "CASSIA",
  "GPTCelltype",
  "mLLMCelltype_claude_sonnet4.5",
  "mLLMCelltype_gpt5",
  "mLLMCelltype_gemini_pro2.5",
  "mLLMCelltype_grok4_fast",
  "mLLMCelltype_deepseek_v3.1_terminus",
  "mLLMCelltype_qwen3_max",
  "mLLMCelltype_llama4_maverick"
)


# --- Marker Based ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = marker_based)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/marker_based_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Marker-based segment: ", e$message)
})

# --- Reference Based ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = reference_based)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/reference_based_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Reference-based segment: ", e$message)
})

# --- DL Based 1 ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = DL_based1)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/DL_based1_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in DL-based segment: ", e$message)
})

# --- Classic ML Based 1 ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = classicML_based1)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/classicML_based1_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Classic ML 1 segment: ", e$message)
})

# --- Classic ML Based 2 ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = classicML_based2)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/classicML_based2_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Classic ML 2 segment: ", e$message)
})

# --- Classic ML Based 3 ---
tryCatch({
  results <- run_all_tools_single_split(seurat_obj, tools_to_run = classicML_based3)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/classicML_based3_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Classic ML 3 segment: ", e$message)
})






# test <- c(
#   "NeuCA"
# )
# 
# results <- run_all_tools_cv(seurat_obj, tools_to_run = test)
# final_results <- aggregate_all_tools(results)
# print(final_results$comparison_table)
# print(final_results$summary_stats)


