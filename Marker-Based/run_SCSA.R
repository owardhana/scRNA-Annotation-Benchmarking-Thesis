# run_SCSA.R
#################################################
# SCSA Function for Benchmarking Framework
# Input: Train/test Seurat objects and markers from FindAllMarkers
# Output: Standardized results format for CV framework
#################################################

#' SCSA Cell Type Annotation Function
#' 
#' Runs SCSA algorithm using conda environment and returns cluster-level predictions
#' 
#' @param seurat_train Training Seurat object (for interface consistency)
#' @param seurat_test Test Seurat object to predict 
#' @param markers Marker genes dataframe from FindAllMarkers()
#' @param scsa_py Path to SCSA.py script
#' @param db_path Path to SCSA database file
#' @param scsa_env Name of conda environment with SCSA
#' @param debug Enable debug output (unused in minimal version)
#' @param keep_files Keep temporary files (default: TRUE)
#' 
#' @return List with predictions, true_labels, confidence_scores, cell_ids
run_SCSA_function <- function(seurat_train, seurat_test, markers,
                              scsa_py = "marker-based/SCSA.py",
                              db_path = "marker-based/whole.db",
                              scsa_env = "scsa_env",
                              debug = FALSE,
                              keep_files = TRUE) {
  
  # Helper function for default return when SCSA fails
  default_return <- function(test_clusters = NULL) {
    if (is.null(test_clusters)) {
      if ("seurat_clusters" %in% colnames(seurat_test@meta.data)) {
        test_clusters <- as.character(seurat_test@meta.data$seurat_clusters)
      } else if ("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
        test_clusters <- as.character(seurat_test@meta.data$Ground_Truth_Celltype)
      } else {
        test_clusters <- rep("Unknown", ncol(seurat_test))
      }
    }
    
    return(list(
      predictions = rep("Unknown", ncol(seurat_test)),
      true_labels = test_clusters,
      confidence_scores = rep(0, ncol(seurat_test)),
      cell_ids = colnames(seurat_test),
      runtime_secs = NA,
      peak_system_memory_mb = NA
    ))
  }
  
  # Validate inputs
  if (!is.data.frame(markers) || nrow(markers) == 0) {
    warning("Invalid markers data")
    return(default_return())
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
  
  # Check conda and environment - prefer miniforge3 installation 
  conda_path <- Sys.which("conda")
  
  # If using homebrew conda, try to use miniforge3 instead for consistency
  if (grepl("homebrew", conda_path)) {
    miniforge_conda <- "/Users/oliverwardhana/miniforge3/bin/conda"
    if (file.exists(miniforge_conda)) {
      conda_path <- miniforge_conda
    }
  }
  
  if (conda_path == "" || !file.exists(conda_path)) {
    warning("conda not found in PATH")
    return(default_return())
  }
  
  # Check if conda environment exists and validate it
  env_check <- tryCatch({
    system2(conda_path, c("env", "list"), stdout = TRUE, stderr = TRUE)
  }, error = function(e) NULL)
  
  env_exists <- FALSE
  if (!is.null(env_check)) {
    env_exists <- any(grepl(paste0("\\b", scsa_env, "\\b"), env_check))
  }
  
  if (!env_exists) {
    warning("Conda environment '", scsa_env, "' not found. Please create it first.")
    return(default_return())
  }
  
  tryCatch({
  # Create temporary directory
  tmpd <- tempfile("SCSA_")
  dir.create(tmpd)
  on.exit({
    if (!keep_files) {
      unlink(tmpd, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  
  # Prepare markers data for SCSA
  markers_to_write <- markers
  
  # Ensure gene column exists
  if (!("gene" %in% tolower(names(markers_to_write)))) {
    if (!is.null(rownames(markers_to_write))) {
      markers_to_write <- cbind(gene = rownames(markers_to_write), markers_to_write)
    } else {
      warning("Could not find gene column in markers")
      return(default_return())
    }
  }
  
  # Fix column names for SCSA compatibility (avg_log2FC -> avg_logFC)
  names(markers_to_write)[tolower(names(markers_to_write)) == "avg_log2fc"] <- "avg_logFC"
  
  # Write markers CSV file
  markers_csv <- file.path(tmpd, "markers.csv")
  utils::write.csv(markers_to_write, file = markers_csv, row.names = FALSE, quote = FALSE)
  
  # Create markers.table file (cluster-gene pairs)
  cluster_col <- NULL
  for (col_name in c("cluster", "ident", "group")) {
    if (col_name %in% tolower(names(markers))) {
      cluster_col <- names(markers)[which(tolower(names(markers)) == col_name)[1]]
      break
    }
  }
  
  cluster_vals <- if (is.null(cluster_col)) rep("cluster1", nrow(markers)) else as.character(markers[[cluster_col]])
  gene_vals <- as.character(markers_to_write$gene)
  markers_table <- data.frame(Cluster = cluster_vals, Gene = gene_vals, stringsAsFactors = FALSE)
  
  markers_table_path <- file.path(tmpd, "markers.table")
  utils::write.table(markers_table, file = markers_table_path, sep = "\t", 
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Execute SCSA via conda
  out_prefix <- file.path(tmpd, "result")
  scsa_args <- c("-d", db_path_abs,
                 "-i", normalizePath(markers_csv),
                 "-s", "seurat",
                 "-o", out_prefix,
                 "-m", "txt",
                 "-M", normalizePath(markers_table_path),
                 "-f", "0.5",    # STANDARDIZED: foldchange threshold (avg_log2FC >= 0.5)
                 "-p", "0.05",   # STANDARDIZED: p-value threshold (p_val_adj < 0.05)
                 "-E", "-N", "-b")  # Gene symbol, no reference DB, no print flags
  
  conda_cmd <- c("run", "-n", scsa_env, "python3", scsa_py_abs, scsa_args)
  
  scsa_output <- tryCatch({
    system2(conda_path, conda_cmd, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    warning("SCSA execution failed: ", e$message)
    return(NULL)
  })

  # Check if SCSA execution succeeded
  exit_status <- attr(scsa_output, "status")
  if (!is.null(exit_status) && exit_status != 0) {
    warning("SCSA failed with exit status: ", exit_status)
    return(default_return())
  }

  # Extract runtime and peak system memory from Python stdout
  runtime_secs <- NA
  peak_system_memory_mb <- NA
  if (!is.null(scsa_output)) {
    runtime_line <- grep("RUNTIME_SECS:", scsa_output, value = TRUE)
    if (length(runtime_line) > 0) {
      runtime_secs <- as.numeric(sub(".*RUNTIME_SECS:([0-9.]+).*", "\\1", runtime_line[1]))
    }
    memory_line <- grep("PEAK_SYSTEM_MEMORY_MB:", scsa_output, value = TRUE)
    if (length(memory_line) > 0) {
      peak_system_memory_mb <- as.numeric(sub(".*PEAK_SYSTEM_MEMORY_MB:([0-9.]+).*", "\\1", memory_line[1]))
    }
  }
  
  # Parse SCSA results
  possible_result_files <- c(paste0(out_prefix, ".txt"), out_prefix, paste0(out_prefix, ".tsv"))
  
  scsa_results <- NULL
  for (rf in possible_result_files) {
    if (file.exists(rf)) {
      scsa_results <- tryCatch({
        utils::read.table(rf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      if (!is.null(scsa_results) && nrow(scsa_results) > 0) break
    }
  }
  
  # Get cluster assignments for return format
  if ("seurat_clusters" %in% colnames(seurat_test@meta.data)) {
    test_clusters <- as.character(seurat_test@meta.data$seurat_clusters)
  } else if ("Ground_Truth_Celltype" %in% colnames(seurat_test@meta.data)) {
    test_clusters <- as.character(seurat_test@meta.data$Ground_Truth_Celltype)
  } else {
    # Fallback: create clustering
    seurat_test <- FindNeighbors(seurat_test, verbose = FALSE)
    seurat_test <- FindClusters(seurat_test, verbose = FALSE)
    test_clusters <- as.character(seurat_test@meta.data$seurat_clusters)
  }
  
  # Check if parsing failed
  if (is.null(scsa_results) || nrow(scsa_results) == 0) {
    warning("Could not parse SCSA results")
    return(default_return(test_clusters))
  }
  
  # Extract top prediction per cluster (highest Z-score)
  scsa_top <- do.call(rbind, lapply(split(scsa_results, scsa_results$Cluster), function(df) {
    df[which.max(df$Z.score), ]
  }))
  
  # Return standardized format for benchmarking framework
  return(list(
    predictions = as.character(scsa_top$Cell.Type),
    true_labels = as.character(scsa_top$Cluster),
    confidence_scores = scsa_top$Z.score,
    cell_ids = colnames(seurat_test),
    runtime_secs = runtime_secs,
    peak_system_memory_mb = peak_system_memory_mb
  ))
  }, error = function(e) {
    warning(sprintf("SCSA failed: %s", e$message))
    return(default_return())
  })
}

# For backward compatibility
run_SCSA <- run_SCSA_function