# ============================================================
# Regenerate the missing singular ground-truth UMAP for
# Nowakowski Cortex 2017 (S3, real arm) — for Supplementary Fig. S4.
# Replicates the UMAP block of P5/Plots_Analysis.R (lines ~2557-2581)
# using the Seurat object from Data Creation (the object was absent
# from benchmarking/data, which is why the singular file was skipped).
# Output filename matches its 8 siblings so it drops straight into S4.
# ============================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

# Both paths point at data artifacts not included in this archive (the preprocessed Seurat
# object from 01-data-generation/real-data-preprocessing/ and the results/ plots output dir).
RDS  <- Sys.getenv("NOWAKOWSKI_RDS", "Nowakowski-Cortex-2017_for_use_subset.rds")
OUT  <- Sys.getenv("NOWAKOWSKI_UMAP_OUT", "Nowakowski-Cortex-2017_for_use_subset_UMAP.png")

seurat_obj <- readRDS(RDS)
stopifnot("Ground_Truth_Celltype" %in% colnames(seurat_obj@meta.data))

Idents(seurat_obj) <- "Ground_Truth_Celltype"
n_pcs <- min(30L, ncol(Embeddings(seurat_obj, "pca")))
seurat_obj <- RunUMAP(seurat_obj, dims = seq_len(n_pcs), verbose = FALSE)

umap_truth <- DimPlot(seurat_obj, reduction = "umap",
                      group.by = "Ground_Truth_Celltype", label = FALSE) +
  ggtitle("Nowakowski Cortex 2017\nGround Truth Cell Types")

ggsave(OUT, plot = umap_truth, width = 178, height = 152, units = "mm", dpi = 1200)
cat("Saved:", OUT, "\n")
cat("cells:", ncol(seurat_obj), "| cell types:", length(unique(seurat_obj$Ground_Truth_Celltype)), "\n")
