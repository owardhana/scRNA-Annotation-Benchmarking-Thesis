library(BiocManager)
library(scRNAseq)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(Matrix)
library(dplyr) 

# ==============================================================================
# 0. Setup Directories
# ==============================================================================
dir.create("./original", showWarnings = FALSE, recursive = TRUE)
dir.create("./for_use", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load and Convert
# ==============================================================================
# 1. Load the data
counts_mtx <- readMM("./counts.mtx")
counts_mtx <- as(counts_mtx, "CsparseMatrix") # Optimize for Seurat

genes <- read.table("./genes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
barcodes <- read.table("./barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1

# 2. Name the matrix
rownames(counts_mtx) <- make.unique(as.character(genes))
colnames(counts_mtx) <- as.character(barcodes)

# 3. Create Object
seurat_obj <- CreateSeuratObject(counts = counts_mtx)
seurat_obj <- UpdateSeuratObject(seurat_obj)

# 4. Add Metadata
metadata_df <- read.table("./metadata.tsv", header = TRUE, sep = "\t", row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_df)

# FIND OUT THE NAME OF THE METADATA COLUMN THAT HAS THE CELL TYPE INFO
head(seurat_obj@meta.data)
cell_type_column <- "celltype"  # Update this if the column name is different

message("Actual cells in original: ", ncol(seurat_obj))

# ==============================================================================
# 2. Extract and Write Original
# ==============================================================================
# Extract the raw counts matrix
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")

# ANNDATA SAFETY PATCH 1: Force strict sparse matrix to prevent 0-byte ghost files
counts_matrix <- as(counts_matrix, "CsparseMatrix")

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix, file = "./original/counts.mtx")

# ANNDATA SAFETY PATCH 2: Create a 2-column dataframe for genes.tsv
genes_df_orig <- data.frame(
  gene_id = rownames(seurat_obj),
  gene_symbol = rownames(seurat_obj)
)

# Write the Genes (Features)
write.table(genes_df_orig, file = "./original/genes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Write the Barcodes (Cells)
write.table(colnames(seurat_obj), file = "./original/barcodes.tsv", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Extract and write metadata cleanly using dplyr
metadata_orig <- data.frame(
  barcode = rownames(seurat_obj@meta.data),
  Ground_Truth_Celltype = seurat_obj@meta.data[[cell_type_column]]  
)

# Write metadata
write.table(metadata_orig, file = "./original/metadata.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# ==============================================================================
# 3. Create the Subset (Dynamic "Zone" Logic)
# ==============================================================================
n_cells <- ncol(seurat_obj)

# Determine the target size based on the dataset's natural scale
if (n_cells < 2000) {
  target <- 500
} else if (n_cells < 9000) {
  target <- 3000
} else {
  target <- 15000
}

# Apply 15% tolerance rule
if (n_cells > (target * 1.15)) {
  message("Dataset size (", n_cells, ") exceeds 15% tolerance for target (", target, "). Subsetting...")
  
  # Calculate the sampling fraction (Target / Current Total)
  sampling_fraction <- target / n_cells
  
  # Perform the stratified sampling
  sampled_barcodes <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("tmp_rn") %>%
    group_by(across(all_of(cell_type_column))) %>%
    slice_sample(prop = sampling_fraction) %>%
    ungroup() %>%
    pull(tmp_rn)
  
  # Create the subset
  seurat_subset <- subset(seurat_obj, cells = sampled_barcodes)
  
} else {
  message("Dataset size (", n_cells, ") is within 15% of target (", target, "). Skipping downsampling.")
  # Pass the original object forward so the pipeline doesn't break
  seurat_subset <- seurat_obj
}

# ==============================================================================
# 4. Extract and Write Subset (or Original if skipped)
# ==============================================================================
# Extract the raw counts matrix from the SUBSET
counts_matrix_subset <- GetAssayData(seurat_subset, layer = "counts")

# ANNDATA SAFETY PATCH 1: Force strict sparse matrix
counts_matrix_subset <- as(counts_matrix_subset, "CsparseMatrix")

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix_subset, file = "./for_use/counts.mtx")

# ANNDATA SAFETY PATCH 2: Create a 2-column dataframe for genes.tsv
genes_df_sub <- data.frame(
  gene_id = rownames(seurat_subset),
  gene_symbol = rownames(seurat_subset)
)

# Write the Genes (Features)
write.table(genes_df_sub, file = "./for_use/genes.tsv", 
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
message("Target threshold determined: ", target)
message("Actual cells in original: ", ncol(seurat_obj))
print(table(seurat_obj[[cell_type_column]]))

message("\nActual cells exported for use: ", ncol(seurat_subset))
print(table(seurat_subset[[cell_type_column]]))