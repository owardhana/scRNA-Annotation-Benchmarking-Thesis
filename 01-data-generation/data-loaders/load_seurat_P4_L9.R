# Load L9 Benchmark Simulations into Seurat
# ============================================================================
#
# PURPOSE:
#   Loads all 27 L9 benchmark datasets (9 scenarios × 3 replicates) into
#   Seurat objects and saves both a subsetted (HVG) and full-gene version
#   of each.
#
# FOLDER STRUCTURE EXPECTED:
#   l9_benchmark_sets/
#     S1/
#       rep1/  ← counts.mtx, genes.txt, cells.txt, metadata.csv
#       rep2/
#       rep3/
#     S2/
#       rep1/
#       ...
#     S9/
#       rep3/
#
# OUTPUT (saved inside each rep folder):
#   l9_benchmark/S1/rep1/S1_rep1_seurat_subset.rds
#   l9_benchmark/S1/rep1/S1_rep1_seurat_full.rds
#
# ============================================================================

library(Seurat)
library(Matrix)

# ============================================================================
# CONFIGURATION
# ============================================================================

base_dir <- "l9_benchmark_sets"

# ============================================================================
# DISCOVER ALL REPLICATE DIRECTORIES
# ============================================================================
#
# The L9 structure is two levels deep: base_dir/scenario/replicate/
# list.dirs with recursive=TRUE finds all subdirectories at any depth.
# We then filter to only those that actually contain counts.mtx — this
# means we get the rep-level directories (not the scenario-level ones).

all_dirs  <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
rep_dirs  <- all_dirs[file.exists(file.path(all_dirs, "counts.mtx"))]

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  L9 Benchmark — Load to Seurat                                    ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("\nBase directory : %s\n", base_dir))
cat(sprintf("Datasets found : %d  (expect 27)\n\n", length(rep_dirs)))

# ============================================================================
# PROCESSING LOOP
# ============================================================================

for (rep_dir in sort(rep_dirs)) {

  # Derive a human-readable name from the path structure:
  # l9_benchmark/S1/rep1 → "S1_rep1"
  parts    <- strsplit(rep_dir, .Platform$file.sep)[[1]]
  sim_name <- paste(tail(parts, 2), collapse = "_")   # e.g. "S1_rep1"

  cat(rep("=", 70), "\n", sep = "")
  cat("Processing:", sim_name, "\n")
  cat("Path:      ", rep_dir, "\n")
  cat(rep("=", 70), "\n", sep = "")



  # --- 1. LOAD DATA -------------------------------------------------------

  counts   <- readMM(file.path(rep_dir, "counts.mtx"))
  genes    <- readLines(file.path(rep_dir, "genes.txt"))
  cells    <- readLines(file.path(rep_dir, "cells.txt"))
  metadata <- read.csv(file.path(rep_dir, "metadata.csv"), row.names = 1)

  rownames(counts) <- genes
  colnames(counts) <- cells

  # Ground truth labels — consistent column name used by benchmarking framework
  metadata$Ground_Truth_Celltype <- metadata$cell_type
  unique_types                   <- sort(unique(metadata$cell_type))
  type_to_cluster                <- setNames(seq_along(unique_types), unique_types)
  metadata$Ground_Truth_Cluster  <- type_to_cluster[metadata$cell_type]



  # --- 2. GLOBAL PREPROCESSING --------------------------------------------

  # Create Seurat object with no automatic filters (applied manually below
  # to match the exact sequential logic used in the scanpy pipeline)
  seurat <- CreateSeuratObject(
    counts    = counts,
    meta.data = metadata,
    project   = sim_name
  )

  # Step A: Filter cells — matches sc.pp.filter_cells(min_genes=200)
  seurat <- subset(seurat, subset = nFeature_RNA >= 200)

  # Step B: Filter genes — matches sc.pp.filter_genes(min_cells=3)
  # Applied after cell filtering to match scanpy's sequential logic
  counts_matrix <- GetAssayData(seurat, layer = "counts")
  genes_to_keep <- rownames(counts_matrix)[rowSums(counts_matrix > 0) >= 3]
  seurat        <- subset(seurat, features = genes_to_keep)

  # Step C: Find variable features (seurat_v3 / vst — matches scanpy flavor)
  seurat <- FindVariableFeatures(
    seurat,
    selection.method = "vst",
    nfeatures        = 2000,
    verbose          = FALSE
  )

  # Step D: Normalize (log1p — matches sc.pp.normalize_total + sc.pp.log1p)
  seurat <- NormalizeData(
    seurat,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE
  )

  cat(sprintf("  Cells after QC : %d\n", ncol(seurat)))
  cat(sprintf("  Genes after QC : %d\n", nrow(seurat)))



  # --- 3. FORK A: SUBSETTED (Standard Benchmark — HVGs only) -------------

  cat("  → Creating Subsetted version (2000 HVGs)...\n")

  seurat_sub <- subset(seurat, features = VariableFeatures(seurat))
  seurat_sub <- ScaleData(seurat_sub, scale.max = 10, verbose = FALSE)
  seurat_sub <- RunPCA(seurat_sub, npcs = 50, verbose = FALSE)

  subset_file <- file.path(rep_dir, paste0(sim_name, "_seurat_subset.rds"))
  saveRDS(seurat_sub, subset_file)
  rm(seurat_sub)



  # --- 4. FORK B: FULL OBJECT (All genes — stress test) ------------------

  cat("  → Creating Full version (all genes, scaled + PCA)...\n")

  seurat_full <- ScaleData(
    seurat,
    features  = rownames(seurat),
    scale.max = 10,
    verbose   = FALSE
  )
  seurat_full <- RunPCA(
    seurat_full,
    features = rownames(seurat_full),
    npcs     = 50,
    verbose  = FALSE
  )

  full_file <- file.path(rep_dir, paste0(sim_name, "_seurat_full.rds"))
  saveRDS(seurat_full, full_file)



  # --- 5. CLEAN UP --------------------------------------------------------

  rm(seurat, seurat_full, counts, metadata, counts_matrix)
  gc(verbose = FALSE)

  cat(sprintf("  ✓ Saved: %s\n\n", sim_name))
}

cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║  ✓ All datasets processed                                         ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("Output: %d Seurat objects (subset + full) across %d datasets\n\n",
            length(rep_dirs) * 2, length(rep_dirs)))
