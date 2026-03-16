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

#################################################
# CONFIGURATION SECTION
#################################################

# Data paths - Datasets (Seurat with Normalized and Find Variable Features run)
SEURAT_PATH <- "./data/imbalanced_subtle_8types_seurat_full.rds"

# Cross-validation configuration
CV_CONFIG <- list(
  k = 3,                              # Number of folds
  group_by = NA,              # For grouped CV, set to NA for stratified
  stratify_by = "Ground_Truth_Celltype", # For stratified CV
  seed = 123                          # Reproducibility
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

DL_based <- c(
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
  results <- run_all_tools_cv(seurat_obj, tools_to_run = marker_based)
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
  results <- run_all_tools_cv(seurat_obj, tools_to_run = reference_based)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/reference_based_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Reference-based segment: ", e$message)
})

# --- DL Based ---
tryCatch({
  results <- run_all_tools_cv(seurat_obj, tools_to_run = DL_based)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/DL_based_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in DL-based segment: ", e$message)
})

# --- Classic ML Based 1 ---
tryCatch({
  results <- run_all_tools_cv(seurat_obj, tools_to_run = classicML_based1)
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
  results <- run_all_tools_cv(seurat_obj, tools_to_run = classicML_based2)
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
  results <- run_all_tools_cv(seurat_obj, tools_to_run = classicML_based3)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(final_results$comparison_table, sprintf("recent_run_data/classicML_based3_tool_comparison_%s.csv", timestamp))
}, error = function(e) {
  message("Error in Classic ML 3 segment: ", e$message)
})






test <- c(
  "CALLR"
)

tryCatch({
  results <- run_all_tools_cv(seurat_obj, tools_to_run = test)
  final_results <- aggregate_all_tools(results)
  print(final_results$comparison_table)
  print(final_results$summary_stats)
}, error = function(e) {
  message("Error in test tools segment: ", e$message)
})


results <- run_all_tools_cv(seurat_obj, tools_to_run = 'SCINA')
final_results <- aggregate_all_tools(results)
print(final_results$comparison_table)
print(final_results$summary_stats)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(CALLR$comparison_table, sprintf("recent_run_data/classicML_based3_tool_comparison_%s.csv", timestamp))

