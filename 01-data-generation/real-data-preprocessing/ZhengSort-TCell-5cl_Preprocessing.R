library(BiocManager)
library(scRNAseq)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr) 

# ==============================================================================
# 0. Setup Directories
# ==============================================================================
dir.create("./original", showWarnings = FALSE, recursive = TRUE)
dir.create("./for_use", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load and Convert
# ==============================================================================
# expr_df: rows = cells, cols = genes
expr_df  <- read.csv("Filtered_DownSampled_SortedPBMC_data.csv", row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)
labels_df <- read.csv("Labels.csv", header = TRUE, stringsAsFactors = FALSE)

# sanity
stopifnot(nrow(expr_df) == nrow(labels_df))

# make counts matrix (genes x cells) for Seurat
counts <- t(as.matrix(expr_df))
cell_names <- rownames(expr_df)
colnames(counts) <- cell_names
rownames(counts) <- colnames(expr_df)

# create metadata data.frame with correct rownames and one column named Ground_Truth_Celltype
meta_df <- data.frame(Ground_Truth_Celltype = as.character(labels_df[[1]]),
                      row.names = cell_names,
                      stringsAsFactors = FALSE)

# create Seurat object
seurat <- CreateSeuratObject(counts = counts, meta.data = meta_df, min.cells = 3,   
                             min.features = 200)

# Make a factor of unique cell types
celltype_factor <- factor(seurat$Ground_Truth_Celltype)

# Convert factor to numeric (1,2,3,...)
seurat$Ground_Truth_Cluster <- as.numeric(celltype_factor)

# Check mapping
levels(celltype_factor)   # shows which celltype maps to which number
head(seurat@meta.data)
seurat_obj <- seurat
seurat_obj <- UpdateSeuratObject(seurat_obj) 
cell_type_column <- "Ground_Truth_Celltype"  # Update this if the column name is different

message("Actual cells in original: ", ncol(seurat_obj))

# ==============================================================================
# 2. Extract and Write Original
# ==============================================================================
# Extract the raw counts matrix
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix, file = "./original/counts.mtx")

# Write the Genes (Features)
write.table(rownames(seurat_obj), file = "./original/genes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Write the Barcodes (Cells)
write.table(colnames(seurat_obj), file = "./original/barcodes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Extract and write metadata cleanly using base R
metadata_orig <- data.frame(
  barcode = rownames(seurat_obj@meta.data),
  Ground_Truth_Celltype = seurat_obj@meta.data[[cell_type_column]]  
)

# Write metadata
write.table(metadata_orig, file = "./original/metadata.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# ==============================================================================
# 3. Create the Subset (Filter for Capital "T")
# ==============================================================================
# Get all unique clusters in the dataset
unique_clusters <- unique(seurat_obj@meta.data[[cell_type_column]])

# Filter for clusters that contain a capital "T"
target_clusters <- grep("T", unique_clusters, value = TRUE)

# Safety check: ensure we actually found some matches!
if(length(target_clusters) == 0) {
  stop("Error: No clusters containing a capital 'T' were found in the dataset.")
}

message("\nSubsetting to keep only the following clusters containing 'T': \n  - ", paste(target_clusters, collapse = "\n  - "))

# Identify all cell barcodes that belong to these specific clusters
cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cell_type_column]] %in% target_clusters]

# Create the subset using those barcodes
seurat_subset <- subset(seurat_obj, cells = cells_to_keep)


# ==============================================================================
# 4. Extract and Write Subset
# ==============================================================================
# Extract the raw counts matrix from the SUBSET
counts_matrix_subset <- GetAssayData(seurat_subset, layer = "counts")

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix_subset, file = "./for_use/counts.mtx")

# Write the Genes (Features)
write.table(rownames(seurat_subset), file = "./for_use/genes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Write the Barcodes (Cells)
write.table(colnames(seurat_subset), file = "./for_use/barcodes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Extract metadata for the subset
metadata_subset <- data.frame(
  barcode = rownames(seurat_subset@meta.data),
  Ground_Truth_Celltype = seurat_subset@meta.data[[cell_type_column]]
)

# Write metadata for the subset
write.table(metadata_subset, file = "./for_use/metadata.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# Verification
# ==============================================================================
message("\n--- Summary ---")
message("Actual cells in original: ", ncol(seurat_obj))
print(table(seurat_obj[[cell_type_column]]))

message("\nActual cells exported for use (Filtered for 'T'): ", ncol(seurat_subset))
print(table(seurat_subset[[cell_type_column]]))