# run_SCSA.R
#################################################
# SCSA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#
# NOTE: SCSA is a database-matching tool, not a supervised classifier.
# It receives marker genes for each cluster and looks them up against a
# curated database to identify cell types. It cannot take reference markers
# and classify individual test cells. Therefore:
#   - seurat_train and markers (training-side) are NOT used.
#   - Test cells are clustered unsupervised, markers are found on the test
#     clusters, SCSA annotates those clusters, and predictions are mapped
#     back to individual cells.
#################################################

#' SCSA Cell Type Annotation Function
#'
#' Runs SCSA algorithm using conda environment and returns cell-level predictions
#'
#' @param seurat_train Training Seurat object (accepted for interface consistency; not used)
#' @param seurat_test Test Seurat object to predict
#' @param markers Marker genes dataframe from FindAllMarkers() (accepted for interface consistency; not used)
#' @param scsa_py Path to SCSA.py script
#' @param db_path Path to SCSA database file
#' @param scsa_env Name of conda environment with SCSA
#' @param keep_files Keep temporary files (default: TRUE)
#'
#' @return List with predictions, true_labels, confidence_scores, cell_ids
run_SCSA_function <- function(seurat_train, seurat_test, markers,
                              scsa_py = "marker-based/SCSA.py",
                              db_path = "marker-based/whole.db",
                              scsa_env = "scsa_env",
                              keep_files = TRUE) {

  # Default return — always cell-level, always against ground truth
  default_return <- function() {
    return(list(
      predictions           = rep("Unknown", ncol(seurat_test)),
      true_labels           = seurat_test$Ground_Truth_Celltype,
      confidence_scores     = rep(0, ncol(seurat_test)),
      cell_ids              = colnames(seurat_test),
      runtime_secs          = NA,
      peak_system_memory_mb = NA
    ))
  }

  # Validate ground truth
  if (!"Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    warning("Ground_Truth_Celltype not found in test metadata")
    return(list(
      predictions           = rep("Unknown", ncol(seurat_test)),
      true_labels           = rep(NA, ncol(seurat_test)),
      confidence_scores     = rep(0, ncol(seurat_test)),
      cell_ids              = colnames(seurat_test),
      runtime_secs          = NA,
      peak_system_memory_mb = NA
    ))
  }

  # Setup file paths
  scsa_py_abs <- normalizePath(scsa_py, mustWork = FALSE)
  db_path_abs <- normalizePath(db_path, mustWork = FALSE)

  if (!file.exists(scsa_py_abs)) {
    warning("SCSA.py not found at: ", scsa_py_abs)
    return(default_return())
  }

  if (!file.exists(db_path_abs)) {
    warning("SCSA database not found at: ", db_path_abs)
    return(default_return())
  }

  # Locate conda — prefer miniforge3 over homebrew conda
  conda_path <- Sys.which("conda")
  if (grepl("homebrew", conda_path)) {
    miniforge_conda <- path.expand("~/miniforge3/bin/conda")
    if (file.exists(miniforge_conda)) conda_path <- miniforge_conda
  }
  if (conda_path == "" || !file.exists(conda_path)) {
    warning("conda not found in PATH")
    return(default_return())
  }

  # Check conda environment
  env_check <- tryCatch(
    system2(conda_path, c("env", "list"), stdout = TRUE, stderr = TRUE),
    error = function(e) NULL
  )
  if (is.null(env_check) || !any(grepl(paste0("\\b", scsa_env, "\\b"), env_check))) {
    warning("Conda environment '", scsa_env, "' not found. Please create it first.")
    return(default_return())
  }

  tryCatch({
    # Create temporary directory
    tmpd <- tempfile("SCSA_")
    dir.create(tmpd)
    on.exit({
      if (!keep_files) unlink(tmpd, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    runtime_secs <- NA
    peak_system_memory_mb <- NA

    # Extract Ground_Truth_Celltype as cluster labels before entering peakRAM scope
    # (avoids <<- scoping issues inside peakRAM::peakRAM)
    cat("SCSA: using Ground_Truth_Celltype as cluster assignments...\n")
    test_clusters <- as.character(seurat_test@meta.data[["Ground_Truth_Celltype"]])

    # Prepare training markers for user reference DB (-M): filter then extract cluster/gene columns
    source("benchmarking_helpers.R", local = TRUE)
    train_markers_filtered <- prepare_markers(markers)
    if (is.null(train_markers_filtered) || nrow(train_markers_filtered) == 0) {
      warning("No training markers available for SCSA user reference DB")
      return(default_return())
    }
    cat(sprintf("SCSA: %d training marker genes across %d cell types for reference DB\n",
                nrow(train_markers_filtered), length(unique(train_markers_filtered$cluster))))
    cat(sprintf("SCSA: %d test clusters across %d cells\n",
                length(unique(test_clusters)), ncol(seurat_test)))

    if (!requireNamespace("peakRAM", quietly = TRUE)) {
      warning("peakRAM package not available for memory/time tracking")
      start_time <- Sys.time()

      # --- SCSA annotation steps (test-side) ---
      markers_to_write <- .scsa_find_test_markers(seurat_test)
      if (is.null(markers_to_write)) return(default_return())
      scsa_top <- .scsa_run(markers_to_write, train_markers_filtered, tmpd, scsa_py_abs, db_path_abs,
                            conda_path, scsa_env)
      # ----------------------------------------

      runtime_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    } else {
      library(peakRAM)
      peakRAM_result <- peakRAM::peakRAM({
        markers_to_write <<- .scsa_find_test_markers(seurat_test)
        scsa_top <<- if (!is.null(markers_to_write)) {
          .scsa_run(markers_to_write, train_markers_filtered, tmpd, scsa_py_abs, db_path_abs,
                    conda_path, scsa_env)
        } else NULL
      })
      runtime_secs <- peakRAM_result$Elapsed_Time_sec[1]
      peak_system_memory_mb <- peakRAM_result$Peak_RAM_Used_MiB[1]
    }

    if (is.null(scsa_top) || nrow(scsa_top) == 0) {
      warning("SCSA returned no results")
      return(default_return())
    }

    # Parse runtime / memory emitted by SCSA.py (overrides peakRAM if available)
    scsa_stdout <- attr(scsa_top, "stdout")
    if (!is.null(scsa_stdout)) {
      rt_line  <- grep("RUNTIME_SECS:", scsa_stdout, value = TRUE)
      mem_line <- grep("PEAK_SYSTEM_MEMORY_MB:", scsa_stdout, value = TRUE)
      if (length(rt_line)  > 0) runtime_secs         <- as.numeric(sub(".*RUNTIME_SECS:([0-9.]+).*",          "\\1", rt_line[1]))
      if (length(mem_line) > 0) peak_system_memory_mb <- as.numeric(sub(".*PEAK_SYSTEM_MEMORY_MB:([0-9.]+).*", "\\1", mem_line[1]))
    }

    # Map cluster-level SCSA predictions back to individual cells
    # scsa_top$Cluster = numeric cluster IDs ("0", "1", ...)
    # scsa_top$Cell.Type = SCSA's predicted cell type label
    cluster_to_type   <- setNames(as.character(scsa_top$Cell.Type), as.character(scsa_top$Cluster))
    cluster_to_zscore <- setNames(as.numeric(scsa_top$Z.score),     as.character(scsa_top$Cluster))

    cell_predictions <- cluster_to_type[test_clusters]
    cell_predictions[is.na(cell_predictions)] <- "Unknown"

    cell_confidence <- cluster_to_zscore[test_clusters]
    cell_confidence[is.na(cell_confidence)] <- 0

    cat(sprintf("SCSA prediction summary: %d cells, %d Unknown, %d unique types\n",
                length(cell_predictions),
                sum(cell_predictions == "Unknown"),
                length(unique(cell_predictions[cell_predictions != "Unknown"]))))
    print(table(cell_predictions))

    return(list(
      predictions           = as.character(cell_predictions),
      true_labels           = seurat_test$Ground_Truth_Celltype,
      confidence_scores     = as.numeric(cell_confidence),
      cell_ids              = colnames(seurat_test),
      runtime_secs          = runtime_secs,
      peak_system_memory_mb = peak_system_memory_mb
    ))

  }, error = function(e) {
    warning(sprintf("SCSA failed: %s", e$message))
    return(default_return())
  })
}

# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

# (Clustering step moved into run_SCSA_function body to avoid peakRAM <<- scoping issues)

# Step 2: Find markers for test clusters; return prepare_markers()-filtered df
.scsa_find_test_markers <- function(seurat_test) {
  source("benchmarking_helpers.R", local = TRUE)
  Idents(seurat_test) <- seurat_test@meta.data[["Ground_Truth_Celltype"]]
  raw_markers <- tryCatch({
    FindAllMarkers(seurat_test, only.pos = TRUE, verbose = FALSE,
                   min.cells.group = 3)
  }, error = function(e) {
    warning("FindAllMarkers on test clusters failed: ", e$message)
    return(data.frame())
  })
  if (nrow(raw_markers) == 0) {
    warning("No markers found for test clusters")
    return(NULL)
  }
  filtered <- prepare_markers(raw_markers)
  if (nrow(filtered) == 0) {
    warning("No markers remain after prepare_markers() filtering")
    return(NULL)
  }
  cat(sprintf("SCSA: %d marker genes across %d test clusters\n",
              nrow(filtered), length(unique(filtered$cluster))))
  return(filtered)
}

# Step 3: Write marker files, execute SCSA, parse and return scsa_top
# train_markers: filtered markers from training data (used as -M user reference DB)
# markers_to_write: markers from test clusters (used as -i query input)
.scsa_run <- function(markers_to_write, train_markers, tmpd, scsa_py_abs, db_path_abs,
                      conda_path, scsa_env) {
  # Fix column name for SCSA compatibility (avg_log2FC -> avg_logFC)
  names(markers_to_write)[tolower(names(markers_to_write)) == "avg_log2fc"] <- "avg_logFC"

  markers_csv <- file.path(tmpd, "markers.csv")
  utils::write.csv(markers_to_write, file = markers_csv, row.names = FALSE, quote = FALSE)

  # markers.table: training markers as user reference DB (cluster-gene pairs, no header)
  # Using training data avoids circular leakage from scoring test markers against themselves
  train_table <- data.frame(
    Cluster = as.character(train_markers$cluster),
    Gene    = as.character(train_markers$gene),
    stringsAsFactors = FALSE
  )
  markers_table_path <- file.path(tmpd, "markers.table")
  utils::write.table(train_table, file = markers_table_path, sep = "\t",
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

  out_prefix <- file.path(tmpd, "result")
  scsa_args  <- c("-d", db_path_abs,
                  "-i", normalizePath(markers_csv),
                  "-s", "seurat",
                  "-o", out_prefix,
                  "-m", "txt",
                  "-M", normalizePath(markers_table_path),
                  "-f", "0",
                  "-E", "-N", "-b")
  conda_cmd  <- c("run", "-n", scsa_env, "python3", scsa_py_abs, scsa_args)

  scsa_output <- tryCatch(
    system2(conda_path, conda_cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) { warning("SCSA execution failed: ", e$message); NULL }
  )

  exit_status <- attr(scsa_output, "status")
  if (!is.null(exit_status) && exit_status != 0) {
    warning("SCSA exited with status: ", exit_status)
    return(NULL)
  }

  # Parse result file
  possible_files <- c(paste0(out_prefix, ".txt"), out_prefix, paste0(out_prefix, ".tsv"))
  scsa_results <- NULL
  for (rf in possible_files) {
    if (file.exists(rf)) {
      scsa_results <- tryCatch(
        utils::read.table(rf, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(scsa_results) && nrow(scsa_results) > 0) break
    }
  }

  if (is.null(scsa_results) || nrow(scsa_results) == 0) {
    warning("Could not parse SCSA results")
    return(NULL)
  }

  # One row per cluster: take highest Z-score prediction
  scsa_top <- do.call(rbind, lapply(split(scsa_results, scsa_results$Cluster), function(df) {
    df[which.max(df$Z.score), ]
  }))

  # Attach stdout for runtime/memory parsing by the caller
  attr(scsa_top, "stdout") <- scsa_output
  return(scsa_top)
}

# For backward compatibility
run_SCSA <- run_SCSA_function
