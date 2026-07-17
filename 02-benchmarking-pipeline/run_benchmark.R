# Benchmarking.R
#################################################
# Scalable Benchmarking Framework for Cell Type Annotation Tools
#
# Structure:
#   1. Library loading
#   2. Source helpers (benchmarking_helpers.R) and tool runners
#   3. Tool registry (maps tool names -> functions)
#   4. Tool group definitions (convenience vectors)
#   5. Configuration (paths, CV settings)
#   6. Example usage for each execution mode
#
# Execution mode functions (defined in benchmarking_helpers.R):
#   - run_all_tools_cv()                  K-fold cross-validation
#   - run_all_tools_single_split()        Single 80/20 train/test split
#   - run_all_tools_full_dataset_llm()    Full dataset (LLM/marker tools)
#   - run_all_tools_full_dataset_database()  Full dataset (CellMarker DB)
#   - run_separate_datasets()             Cross-dataset generalization
#
# All outputs are compatible with aggregate_all_tools() for comparison.
#################################################
# Load required libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(caret)
library(SingleCellExperiment)
library(SummarizedExperiment)
# library(monocle)
library(tensorflow)
library(reticulate)
library(bench)
library(pryr)
suppressWarnings(library(logr))

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
source("marker-based/run_CellAssign.R")

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

source("database-based/run_SCINA.R")
source("database-based/run_scSorter.R")
source("database-based/run_scType.R")
source("database-based/run_scCATCH.R")
source("database-based/run_Garnett.R")
source("database-based/run_SCSA.R")
source("database-based/run_clustifyr_hyper.R")
source("database-based/run_clustifyr_jaccard.R")
source("database-based/run_ScInfeR.R")


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
TOOL_REGISTRY[["CellAssign"]] <- run_CellAssign_function

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
TOOL_REGISTRY[["scPred_svmLinear"]] <- run_scPred_svmLinear_function
TOOL_REGISTRY[["scPred_nb"]] <- run_scPred_nb_function
TOOL_REGISTRY[["scPred_svmPoly"]] <- run_scPred_svmPoly_function
TOOL_REGISTRY[["scPred_bayesglm"]] <- run_scPred_bayesglm_function
TOOL_REGISTRY[["scPred_earth"]] <- run_scPred_earth_function
TOOL_REGISTRY[["scPred_mlp"]] <- run_scPred_mlp_function
TOOL_REGISTRY[["scPred_nnet"]] <- run_scPred_nnet_function
TOOL_REGISTRY[["scPred_regLogistic"]] <- run_scPred_regLogistic_function
TOOL_REGISTRY[["scPred_multinom"]] <- run_scPred_multinom_function
TOOL_REGISTRY[["scPred_glmnet"]] <- run_scPred_glmnet_function

TOOL_REGISTRY[["scLearn"]] <- run_scLearn_function
TOOL_REGISTRY[["CAMLU"]] <- run_CAMLU_function
TOOL_REGISTRY[["NeuCA"]] <- run_NeuCA_function
TOOL_REGISTRY[["scPred_mxnetAdam"]] <- run_scPred_mxnetAdam_function

TOOL_REGISTRY[["mLLMCelltype"]] <- run_mLLMCelltype_function
TOOL_REGISTRY[["CASSIA"]] <- run_CASSIA_function
TOOL_REGISTRY[["GPTCelltype"]] <- run_GPTCelltype_function

TOOL_REGISTRY[["SCINA_database"]] <- run_SCINA_database_function
TOOL_REGISTRY[["scSorter_database"]] <- run_scSorter_database_function
TOOL_REGISTRY[["scType_database"]] <- run_scType_database_function
TOOL_REGISTRY[["scCATCH_database"]] <- run_scCATCH_database_function
TOOL_REGISTRY[["Garnett_database"]] <- run_Garnett_database_function
TOOL_REGISTRY[["SCSA_database"]] <- run_SCSA_database_function
TOOL_REGISTRY[["clustifyr_hyper_database"]] <- run_clustifyr_hyper_database_function
TOOL_REGISTRY[["clustifyr_jaccard_database"]] <- run_clustifyr_jaccard_database_function
TOOL_REGISTRY[["ScInfeR_database"]] <- run_ScInfeR_database_function




#################################################
# TOOL GROUP DEFINITIONS
#################################################
# Organized by method category for convenient batch runs.
# Uncomment/comment individual tools as needed.

marker_based <- c(
  "SCINA", "scSorter", "scType", "scCATCH", "SCSA",
  # "Garnett",        # broken
  "clustifyr_hyper", "clustifyr_jaccard", "ScInfeR"
  # "CellAssign"      # broken
)

reference_based <- c(
  "SingleR", "Seurat_Transfer_PCA", "Seurat_Transfer_RPCA",
  "Seurat_Transfer_CCA", "scmap_cell", "scmap_cluster",
  "scibetR", "CIPR"
)

classicML_based1 <- c(
  "scPred", "scPred_avNNet", "scPred_xgbTree", "scPred_rf",
  "scPred_glm", "scPred_glmboost", "scPred_lda", "scPred_knn",
  "CHETAH", "scClassify", "scAnnotatR", "singleCellNet", "CellTypist"
)

classicML_based2 <- c(
  "scAnnotate", "CALLR", "CaSTLe",
  # "scPred_adaboost",  # too slow
  "scID"
)

classicML_based3 <- c(
  "scPred_svmLinear", "scPred_nb", "scPred_svmPoly",
  "scPred_bayesglm", "scPred_earth", "scPred_mlp",
  "scPred_nnet", "scPred_regLogistic", "scPred_multinom"
  # "scPred_glmnet"     # broken
)

DL_based <- c(
  "scLearn", "CAMLU", "NeuCA"
  # "scPred_mxnetAdam"  # broken
)

LLM_based <- c(
  "CASSIA", "GPTCelltype", "mLLMCelltype"
)

DB_Tools <- c(
  "SCINA_database", "scType_database", #"Garnett_database",
  "scSorter_database", "scCATCH_database", "SCSA_database",
  "clustifyr_hyper_database", "clustifyr_jaccard_database",
  "ScInfeR_database"
)


#################################################
# EXAMPLE USAGE
#################################################
# Each execution mode below is self-contained. Uncomment the mode you want.
# All outputs are compatible with aggregate_all_tools() for comparison tables.


# -------------------------------------------------
# MODE 1: K-Fold Cross-Validation  (run_all_tools_cv)
#   Robust within-dataset evaluation with fold variation stats.
#   Compatible: marker-based, reference-based, classic-ML, DL
# -------------------------------------------------
# seurat_obj <- readRDS(SEURAT_PATH)
# Idents(seurat_obj) <- "Ground_Truth_Celltype"
#
# CV_CONFIG <- list(
#   k = 3,
#   group_by = NA,                       # Set to "donor_id" for grouped CV
#   stratify_by = "Ground_Truth_Celltype",
#   seed = 123
# )
#
# results <- run_all_tools_cv(seurat_obj,
#                              tools_to_run = reference_based,
#                              cv_config = CV_CONFIG)
# final <- aggregate_all_tools(results)
# print(final$comparison_table)
# write.csv(final$comparison_table,
#           sprintf("recent_run_data/cv_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# -------------------------------------------------
# MODE 2: Single 80/20 Split  (run_all_tools_single_split)
#   Fast single-run evaluation. Good for simulated (L9) datasets.
#   Compatible: marker-based, reference-based, classic-ML, DL
# -------------------------------------------------
# seurat_obj <- readRDS(SEURAT_PATH)
# Idents(seurat_obj) <- "Ground_Truth_Celltype"
#
# SPLIT_CONFIG <- list(
#   label_col = "Ground_Truth_Celltype",
#   train_prop = 0.80,
#   seed = 123
# )
# results <- run_all_tools_single_split(seurat_obj,
#                                        tools_to_run = marker_based)
# final <- aggregate_all_tools(results)
# print(final$comparison_table)
# write.csv(final$comparison_table,
#           sprintf("recent_run_data/single_split_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# -------------------------------------------------
# MODE 3: Full Dataset - LLM/Marker  (run_all_tools_full_dataset_llm)
#   LLM/marker tools with no training step. Uses Ground_Truth_Cluster
#   for perfect clustering + FindAllMarkers on full dataset.
#   Compatible: LLM-based, marker-based (full dataset evaluation)
# -------------------------------------------------
# seurat_obj <- readRDS(SEURAT_PATH)
# Idents(seurat_obj) <- "Ground_Truth_Celltype"
#
# results <- run_all_tools_full_dataset_llm(seurat_obj,
#                                            tools_to_run = LLM_based,
#                                            dataset_name = "Darmanis-Brain-2015")
# final <- aggregate_all_tools(results)
# print(final$comparison_table)
# write.csv(final$comparison_table,
#           sprintf("recent_run_data/llm_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# -------------------------------------------------
# MODE 4: Full Dataset - Database  (run_all_tools_full_dataset_database)
#   Marker-database tools using CellMarker DB instead of FindAllMarkers.
#   Requires tissue_class for DB filtering.
#   Compatible: database-based (*_database tool variants)
# -------------------------------------------------
# seurat_obj <- readRDS(SEURAT_PATH)
# Idents(seurat_obj) <- "Ground_Truth_Celltype"
#
# results <- run_all_tools_full_dataset_database(seurat_obj,
#                                                 tools_to_run  = DB_Tools,
#                                                 tissue_class  = "Brain",
#                                                 dataset_name  = "Darmanis-Brain-2015")
# final <- aggregate_all_tools(results)
# print(final$comparison_table)
# write.csv(final$comparison_table,
#           sprintf("recent_run_data/database_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# -------------------------------------------------
# MODE 5: Separate Datasets  (run_separate_datasets)
#   Train on Dataset A, test on Dataset B. Evaluates cross-dataset
#   generalization - the most realistic annotation use case.
#   Compatible: reference-based, classic-ML, DL
# -------------------------------------------------
# seurat_train <- readRDS("data/Darmanis-Brain-2015_for_use_subset.rds")
# seurat_test  <- readRDS("data/Marques-Brain-2016_for_use_subset.rds")
#
# results <- run_separate_datasets(seurat_train, seurat_test,
#                                   tools_to_run = reference_based)
# final <- aggregate_all_tools(results)
# print(final$comparison_table)
# write.csv(final$comparison_table,
#           sprintf("recent_run_data/cross_dataset_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# -------------------------------------------------
# BATCH EXAMPLE: Loop over multiple datasets (Single Split)
# -------------------------------------------------
# datasets <- c(
#   "Darmanis-Brain-2015_for_use_subset.rds",
#   "He-Skin-2020_for_use_subset.rds",
#   "Grun-Pancreas-2016_for_use_subset.rds"
# )
#
# for (rds_file in datasets) {
#   dataset_name <- tools::file_path_sans_ext(rds_file)
#   cat(sprintf("\n\n====== Dataset: %s ======\n", dataset_name))
#   tryCatch({
#     fold_cache <- list()
#     seurat_obj <- readRDS(file.path("data", rds_file))
#     Idents(seurat_obj) <- "Ground_Truth_Celltype"
#     results <- run_all_tools_single_split(seurat_obj, tools_to_run = reference_based)
#     final <- aggregate_all_tools(results)
#     print(final$comparison_table)
#     timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#     write.csv(final$comparison_table,
#               sprintf("recent_run_data/%s_%s.csv", dataset_name, timestamp))
#   }, error = function(e) {
#     message(sprintf("Error for %s: %s", dataset_name, e$message))
#   })
# }


#################################################
# BATCH CROSS-DATASET TRIALS (User Request) Phase 3
#################################################
# Run all tools from L212-L247 on 9 specified cross-dataset pairings




pairings <- list(
  list(name = "CB-1", train = "Cellbench 10xv2_original_subset.rds", test = "Cellbench Cel-Seq2_original_subset.rds"),
  list(name = "CB-2", train = "Cellbench 10xv2_original_subset.rds", test = "Cellbench Drop-seq_original_subset.rds"),
  list(name = "CB-3", train = "Cellbench Cel-Seq2_original_subset.rds", test = "Cellbench 10xv2_original_subset.rds"),
  list(name = "CB-4", train = "Cellbench Drop-seq_original_subset.rds", test = "Cellbench 10xv2_original_subset.rds"),
  list(name = "PX-1", train = "Baron-Pancreas_original_subset.rds", test = "Segerstolpe-Pancreas-2016_original_subset.rds"),
  list(name = "PX-2", train = "Segerstolpe-Pancreas-2016_original_subset.rds", test = "Baron-Pancreas_original_subset.rds"),
  list(name = "PBMC-A", train = "PBMCBench 10xv2-8cl_original_subset.rds", test = "PBMCBench Drop-seq-8cl_original_subset.rds"),
  list(name = "PBMC-B", train = "PBMCBench 10xv2-8cl_original_subset.rds", test = "PBMCBench InDrops-8cl_original_subset.rds"),
  list(name = "PBMC-C", train = "PBMCBench 10xv2-6cl_original_subset.rds", test = "PBMCBench Seq-well-6cl_original_subset.rds")
)

for (pair in pairings) {
  cat(sprintf("\n\n====== Cross-Dataset Trial: %s ======\n", pair$name))
  tryCatch(
    {
      seurat_train <- readRDS(file.path("data", pair$train))
      seurat_test <- readRDS(file.path("data", pair$test))

      results <- run_separate_datasets(seurat_train, seurat_test, tools_to_run = tools_for_trials)
      final <- aggregate_all_tools(results)

      print(final$comparison_table)

      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      write.csv(
        final$comparison_table,
        sprintf("recent_run_data/cross_dataset_trial_%s_%s.csv", pair$name, timestamp)
      )

      # Clear memory to prevent OOM across loops
      rm(seurat_train, seurat_test, results, final)
      gc()
    },
    error = function(e) {
      message(sprintf("Error for pairing %s: %s", pair$name, e$message))
    }
  )
}

#################################################
# ADDITIONAL TRIALS
#################################################
tools_for_trials <- c(
  marker_based,
  reference_based,
  classicML_based1,
  classicML_based2,
  classicML_based3,
  DL_based
)

seurat_obj <- readRDS("./data/S3_rep2_seurat_subset.rds")

# Remove cell types too rare for cross-validation
min_cells <- 3

ct_counts   <- table(seurat_obj$Ground_Truth_Celltype)
valid_types <- names(ct_counts[ct_counts >= min_cells])
removed     <- names(ct_counts[ct_counts < min_cells])

cat(sprintf("Removing: %s\n", paste(removed, collapse = ", ")))
valid_cells <- colnames(seurat_obj)[seurat_obj$Ground_Truth_Celltype %in% valid_types]
seurat_obj  <- subset(seurat_obj, cells = valid_cells)


seurat_obj <- subset(x = seurat_obj, subset = Ground_Truth_Celltype != "")

valid_cells  <- colnames(seurat_obj)[!is.na(seurat_obj$Ground_Truth_Celltype)]
seurat_obj <- subset(seurat_obj, cells = valid_cells)
Idents(seurat_obj) <- "Ground_Truth_Celltype"

results <- run_all_tools_single_split(seurat_obj,
                                       tools_to_run = tools_for_trials, seed = 124)
final <- aggregate_all_tools(results)
print(final$comparison_table)
write.csv(final$comparison_table,
          sprintf("recent_run_data/single_split_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))


# --- 2. MODE 4: Database with DB_Tools ---
mode4_configs <- list(
  list(file = "Grun-Pancreas-2016_original_subset.rds", tissue_class = c("Pancreas")),
  list(file = "Marques-Brain-2016_for_use_subset.rds", tissue_class = c("Brain")),
  list(file = "Tabula Muris-FACS-3k_for_use_subset.rds", tissue_class = c(
    "Adipose tissue",
    "Bladder",
    "Blood",
    "Bone marrow",
    "Brain",
    "Bronchus",
    "Colon",
    "Epidermis",
    "Heart",
    "Pancreas",
    "Skeletal muscle",
    "Skin",
    "Thymus"
  ))
)

for (cfg in scina_rerun_configs) {
  dataset_name <- tools::file_path_sans_ext(cfg$file)
  cat(sprintf("\n\n====== MODE 4 Dataset: %s ======\n", dataset_name))
  tryCatch({
    seurat_obj <- readRDS(file.path("data", cfg$file))
    Idents(seurat_obj) <- "Ground_Truth_Celltype"
    results <- run_all_tools_full_dataset_database(seurat_obj,
                                                   tools_to_run  = c("SCINA_database"),
                                                   tissue_class  = cfg$tissue_class,
                                                   dataset_name  = dataset_name)
    final <- aggregate_all_tools(results)
    print(final$comparison_table)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(final$comparison_table,
              sprintf("recent_run_data/database_%s_%s.csv", dataset_name, timestamp))
    rm(seurat_obj, results, final)
    gc()
  }, error = function(e) {
    message(sprintf("Error for %s: %s", dataset_name, e$message))
  })
}

# =====================================================================
# MODE 4b: SCINA-ONLY RERUN (Phase 4 / P7a) — after the rm_overlap fix
# ---------------------------------------------------------------------
# Regenerates ONLY the SCINA database results across all 9 P7a datasets.
# The previous P7a run used rm_overlap = TRUE, which empties the heavily
# overlapping CellMarker signatures and aborts SCINA ("'a' must have
# dims > 0") -> all-Unknown on 8/9 datasets. run_SCINA.R now uses
# rm_overlap = FALSE (database-based only). This block re-runs SCINA and
# leaves the other 7 tools' results untouched.
#
# Tissue subsets below match what the other DB tools used in P7a
# (reverse-engineered from their predicted cell-type vocabularies).
#
# HOW TO RUN: select this whole block and run it (or source the script
# with the other MODE blocks commented out). Outputs are written to
#   recent_run_data/<dataset>_SCINA_database_cluster_map.csv
#
# AFTER RUNNING: copy those 9 CSVs over the old ones in
#   results/P7a Marker - Real Database/raw_data/
# then re-run score_cells.py (rescoring) and Plots_Analysis.R.
# =====================================================================
scina_rerun_configs <- list(
  list(file = "Darmanis-Brain-2015_for_use_subset.rds",                tissue_class = c("Brain")),
  list(file = "Marques-Brain-2016_for_use_subset.rds",                 tissue_class = c("Brain")),
  list(file = "Nowakowski-Cortex-2017_for_use_subset.rds",             tissue_class = c("Brain", "Fetal brain")),
  list(file = "Grun-Pancreas-2016_original_subset.rds",                tissue_class = c("Pancreas")),
  list(file = "He-Skin-2020_for_use_subset.rds",                       tissue_class = c("Skin")),
  list(file = "MacParland-Liver-Broad-2018_original_seurat_subset.rds", tissue_class = c("Liver")),
  list(file = "Zhao-Immune-Fine-2020_for_use_subset.rds",              tissue_class = c("Blood")),
  list(file = "Zheng-ZhengSort-5cl-2017_original_subset.rds",          tissue_class = c("Blood")),
  list(file = "Tabula Muris-FACS-3k_for_use_subset.rds",               tissue_class = c(
    "Adipose tissue", "Bladder", "Blood", "Bone marrow", "Brain", "Bronchus",
    "Colon", "Epidermis", "Heart", "Pancreas", "Skeletal muscle", "Skin", "Thymus"
  ))
)

for (cfg in scina_rerun_configs) {
  dataset_name <- tools::file_path_sans_ext(cfg$file)
  cat(sprintf("\n\n====== SCINA RERUN Dataset: %s ======\n", dataset_name))
  tryCatch({
    seurat_obj <- readRDS(file.path("data", cfg$file))
    Idents(seurat_obj) <- "Ground_Truth_Celltype"
    run_all_tools_full_dataset_database(seurat_obj,
                                        tools_to_run  = c("SCINA_database"),
                                        tissue_class  = cfg$tissue_class,
                                        dataset_name  = dataset_name)
    rm(seurat_obj)
    gc()
  }, error = function(e) {
    message(sprintf("Error for %s: %s", dataset_name, e$message))
  })
}

# --- 3. MODE 3: LLM/Marker with LLM_based ---
datasets_mode3 <- c(
  "Zhao-Immune-Fine-2020_for_use_subset.rds"
)

for (rds_file in datasets_mode3) {
  dataset_name <- tools::file_path_sans_ext(rds_file)
  cat(sprintf("\n\n====== MODE 3 Dataset: %s ======\n", dataset_name))
  tryCatch({
    seurat_obj <- readRDS(file.path("data", rds_file))
    Idents(seurat_obj) <- "Ground_Truth_Celltype"
    results <- run_all_tools_full_dataset_llm(seurat_obj,
                                              tools_to_run = LLM_based,
                                              dataset_name = dataset_name)
    final <- aggregate_all_tools(results)
    print(final$comparison_table)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(final$comparison_table,
              sprintf("recent_run_data/llm_%s_%s.csv", dataset_name, timestamp))
    rm(seurat_obj, results, final)
    gc()
  }, error = function(e) {
    message(sprintf("Error for %s: %s", dataset_name, e$message))
  })
}
