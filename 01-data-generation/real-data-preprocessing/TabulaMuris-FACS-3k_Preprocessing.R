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
counts <- ReadMtx(
  mtx            = "counts.mtx",
  cells          = "barcodes.csv",
  features       = "genes.csv",
  cell.column    = 1,
  feature.column = 1,
  skip.cell      = 1,
  skip.feature   = 1,
  mtx.transpose  = TRUE
)

# MAKE GENE NAMES ALL CAPS
rownames(counts) <- toupper(rownames(counts))
seurat <- CreateSeuratObject(counts = counts, project = "TabulaMurisFACS")

meta <- read.csv("metadata.csv", row.names = 1)
shared_cells <- intersect(colnames(seurat), rownames(meta))
seurat <- seurat[, shared_cells]
seurat <- AddMetaData(seurat, metadata = meta[shared_cells, ])

#Remove Rik Genes
seurat_obj <- seurat[!grepl("Rik$", rownames(seurat), ignore.case = TRUE), ]

# Check original vs new number of genes
nrow(seurat)
nrow(seurat_obj)

# FIND OUT THE NAME OF THE METADATA COLUMN THAT HAS THE CELL TYPE INFO
head(seurat_obj@meta.data)
cell_type_column <- "cell_ontology_class"  # Update this if the column name is different

message("Actual cells in original: ", ncol(seurat_obj))

# ==============================================================================
# 2. Extract and Write Original
# ==============================================================================
# Extract the raw counts matrix
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix, file = "./original/counts.mtx")

# Write the Genes (Features) - 2 columns for downstream compatibility
genes_df_orig <- data.frame(
  gene_id = rownames(seurat_obj),
  gene_symbol = rownames(seurat_obj)
)
write.table(genes_df_orig, file = "./original/genes.tsv", 
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
# 3. Create the Subset (Forced to 3k)
# ==============================================================================
n_cells <- ncol(seurat_obj)

# Hardcode the target to 3000, overriding the dynamic sizing logic
target <- 3000

# Apply 15% tolerance rule (will only skip if dataset is < 3450 cells)
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

# Write the Matrix Market file (.mtx)
writeMM(counts_matrix_subset, file = "./for_use/counts.mtx")

# Write the Genes (Features) - 2 columns for downstream compatibility
genes_df_sub <- data.frame(
  gene_id = rownames(seurat_subset),
  gene_symbol = rownames(seurat_subset)
)
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
message("Target threshold forced to: ", target)
message("Actual cells in original: ", ncol(seurat_obj))
print(table(seurat_obj[[cell_type_column]]))

message("\nActual cells exported for use: ", ncol(seurat_subset))
print(table(seurat_subset[[cell_type_column]]))