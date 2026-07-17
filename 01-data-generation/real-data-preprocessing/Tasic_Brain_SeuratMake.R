# ============================================================================
# SINGLE DATASET PROCESSING
# ============================================================================
library(Seurat)
library(Matrix)

# 1. Define your locations based on your current folder
study_dir   <- getwd()                       # e.g., ".../MacParland-Liver-2018"
study_name  <- basename(study_dir)           # e.g., "MacParland-Liver-2018"

# 2. Toggle which subfolder you want to process here
folder_type <- "for_use"                     # Change to "original" if needed
target_dir  <- file.path(study_dir, folder_type)

sim_name    <- paste(study_name, folder_type, sep = "_") 

cat(rep("=", 70), "\n", sep = "")
cat("Processing:", sim_name, "\n")
cat("Target:    ", target_dir, "\n")
cat(rep("=", 70), "\n", sep = "")

# ------------------------------------------------------------------------
# "RESUME" CHECK: Run only if not already processed
# ------------------------------------------------------------------------
subset_file <- file.path(target_dir, paste0(sim_name, "_seurat_subset.rds"))

if (file.exists(subset_file)) {
  cat(sprintf("  -> %s already exists. Skipping processing.\n", basename(subset_file)))
  
} else {
  
  # --- 1. LOAD DATA -------------------------------------------------------
  counts_file <- file.path(target_dir, "counts.mtx")
  if (!file.exists(counts_file)) {
    stop(sprintf("Missing counts.mtx in %s. Check your folder_type toggle.", target_dir))
  }
  
  counts <- readMM(counts_file)
  
  # Force the matrix into Seurat's preferred format BEFORE assigning names
  counts <- as(counts, "CsparseMatrix")
  
  # Load the standard biological outputs (.tsv)
  if (file.exists(file.path(target_dir, "genes.tsv"))) {
    genes    <- read.table(file.path(target_dir, "genes.tsv"), sep="\t", header=FALSE, stringsAsFactors=FALSE)$V1
    cells    <- read.table(file.path(target_dir, "barcodes.tsv"), sep="\t", header=FALSE, stringsAsFactors=FALSE)$V1
    metadata <- read.table(file.path(target_dir, "metadata.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # THE UNDERSCORE FIX: Force standard Seurat formatting to prevent "No Cell Overlap"
    clean_barcodes <- gsub("_", "-", as.character(cells))
    clean_barcodes <- make.unique(clean_barcodes)
    
    if("barcode" %in% colnames(metadata)) {
      metadata$barcode <- clean_barcodes
    }
    rownames(metadata) <- clean_barcodes
    colnames(counts)   <- clean_barcodes
    
  } else {
    # Fallback for raw Splatter text files (just in case)
    genes    <- readLines(file.path(target_dir, "genes.txt"))
    cells    <- readLines(file.path(target_dir, "cells.txt"))
    metadata <- read.csv(file.path(target_dir, "metadata.csv"), row.names = 1)
    colnames(counts) <- as.character(cells)
  }
  
  rownames(counts) <- make.unique(as.character(genes))
  
  # Standardize Ground Truth Metadata
  if (!"Ground_Truth_Celltype" %in% colnames(metadata)) {
    if ("cell_type" %in% colnames(metadata)) {
      metadata$Ground_Truth_Celltype <- metadata$cell_type
    } else if ("cell_ontology_class" %in% colnames(metadata)) {
      metadata$Ground_Truth_Celltype <- metadata$cell_ontology_class
    } else {
      stop("Could not find a valid cell type column in metadata.")
    }
  }
  
  unique_types                   <- sort(unique(metadata$Ground_Truth_Celltype))
  type_to_cluster                <- setNames(seq_along(unique_types), unique_types)
  metadata$Ground_Truth_Cluster  <- type_to_cluster[metadata$Ground_Truth_Celltype]
  
  # --- 2. GLOBAL PREPROCESSING --------------------------------------------
  seurat <- CreateSeuratObject(
    counts    = counts,
    meta.data = metadata,
    project   = sim_name
  )
  
  #seurat <- subset(seurat, subset = nFeature_RNA >= 200)
  
  #counts_matrix <- GetAssayData(seurat, layer = "counts")
  #genes_to_keep <- rownames(counts_matrix)[rowSums(counts_matrix > 0) >= 3]
  #seurat        <- subset(seurat, features = genes_to_keep)
  # Purge unlabeled cells (Matches the Python script)
  seurat <- seurat[, !is.na(seurat$Ground_Truth_Celltype)]
  
  
  seurat <- FindVariableFeatures(
    seurat,
    selection.method = "dispersion",
    nfeatures        = 2000,
    verbose          = FALSE
  )
  
  seurat <- NormalizeData(
    seurat,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE
  )
  
  cat(sprintf("  Cells after QC : %d\n", ncol(seurat)))
  cat(sprintf("  Genes after QC : %d\n", nrow(seurat)))
  
  # --- 3. FORK A: SUBSETTED (Standard Benchmark — HVGs only) -------------
  cat("  -> Creating Subsetted version (2000 HVGs)...\n")
  
  seurat_sub <- subset(seurat, features = VariableFeatures(seurat))
  seurat_sub <- ScaleData(seurat_sub, scale.max = 10, verbose = FALSE)
  
  # Determine a safe number of PCs based on dataset size to avoid irlba SVD errors
  safe_npcs <- min(50, ncol(seurat_sub) - 1)
  seurat_sub <- RunPCA(seurat_sub, npcs = safe_npcs, verbose = FALSE, approx = FALSE)
  
  saveRDS(seurat_sub, subset_file)
  
  # --- 4. CLEAN UP --------------------------------------------------------
  rm(seurat, seurat_sub, counts, metadata)
  gc(verbose = FALSE)
  
  cat(sprintf("  ✓ Saved: %s\n\n", basename(subset_file)))
}

cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  ✓ Single dataset processing complete                              ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
