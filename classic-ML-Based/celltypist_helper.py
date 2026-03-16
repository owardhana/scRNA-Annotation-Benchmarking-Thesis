#!/usr/bin/env python3
"""
CellTypist Helper Script - Implementation
Standalone script for CellTypist training and prediction using file-based communication.

Usage:
python celltypist_helper.py <train_file> <test_file> <labels_file> <genes_file> <output_file> <status_file>

Arguments:
- train_file: CSV file with training expression data (cells x genes)
- test_file: CSV file with test expression data (cells x genes)  
- labels_file: Text file with cell type labels (one per line)
- genes_file: Text file with gene names (one per line)
- output_file: CSV file to write predictions
- status_file: Text file to write execution status
"""

import sys
import os
import pandas as pd
import numpy as np
import warnings
import tracemalloc
warnings.filterwarnings('ignore')

# Import with error handling
try:
    import scanpy as sc
    print("✓ scanpy imported successfully")
except ImportError as e:
    print(f"Error importing scanpy: {e}")
    sys.exit(1)

try:
    import celltypist
    from celltypist import models
    print("✓ celltypist imported successfully")
except ImportError as e:
    print(f"Error importing celltypist: {e}")
    sys.exit(1)

def main():
    """Main execution function"""
    
    # Check command line arguments
    if len(sys.argv) != 7:
        print("Error: Incorrect number of arguments")
        print("Usage: python celltypist_helper.py <train_file> <test_file> <labels_file> <genes_file> <output_file> <status_file>")
        sys.exit(1)
    
    # Parse arguments
    train_file = sys.argv[1]
    test_file = sys.argv[2]
    labels_file = sys.argv[3]
    genes_file = sys.argv[4]
    output_file = sys.argv[5]
    status_file = sys.argv[6]
    
    try:
        # Write initial status
        write_status(status_file, "STARTING")
        
        print("CellTypist Helper - Implementation")
        print(f"Training data: {train_file}")
        print(f"Test data: {test_file}")
        print(f"Labels file: {labels_file}")
        print(f"Genes file: {genes_file}")
        print(f"Output file: {output_file}")
        print(f"Status file: {status_file}")
        
        # Validate input files exist
        input_files = [train_file, test_file, labels_file, genes_file]
        for file_path in input_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")
        
        print("✓ All input files found")
        
        # Load data
        print("Loading data...")
        
        # Load expression matrices
        train_df = pd.read_csv(train_file, index_col=0)
        test_df = pd.read_csv(test_file, index_col=0)
        
        print(f"Training data shape: {train_df.shape}")
        print(f"Test data shape: {test_df.shape}")
        
        # Load labels and genes
        with open(labels_file, 'r') as f:
            labels = [line.strip() for line in f.readlines()]
        
        with open(genes_file, 'r') as f:
            genes = [line.strip() for line in f.readlines()]
        
        print(f"Number of labels: {len(labels)}")
        print(f"Number of genes: {len(genes)}")
        
        # Validate dimensions
        if train_df.shape[0] != len(labels):
            raise ValueError(f"Training data rows ({train_df.shape[0]}) != labels count ({len(labels)})")
        
        if train_df.shape[1] != len(genes):
            raise ValueError(f"Training data columns ({train_df.shape[1]}) != genes count ({len(genes)})")
        
        if test_df.shape[1] != len(genes):
            raise ValueError(f"Test data columns ({test_df.shape[1]}) != genes count ({len(genes)})")
        
        print("✓ Data dimensions validated")
        
        # Create AnnData objects
        print("Creating AnnData objects...")
        
        # Training data
        adata_train = sc.AnnData(
            X=train_df.values.astype(np.float32),
            obs=pd.DataFrame(index=train_df.index),
            var=pd.DataFrame(index=genes)
        )
        adata_train.obs['cell_type'] = labels
        adata_train.obs_names = train_df.index
        adata_train.var_names = genes
        
        # Test data  
        adata_test = sc.AnnData(
            X=test_df.values.astype(np.float32),
            obs=pd.DataFrame(index=test_df.index),
            var=pd.DataFrame(index=genes)
        )
        adata_test.obs_names = test_df.index
        adata_test.var_names = genes
        
        print(f"Training AnnData: {adata_train.shape}")
        print(f"Test AnnData: {adata_test.shape}")
        print("✓ AnnData objects created")
        
        # Normalize data (CellTypist expects log-normalized data)
        print("Normalizing data...")
        
        # Normalize training data
        sc.pp.normalize_total(adata_train, target_sum=1e4)
        sc.pp.log1p(adata_train)
        
        # Normalize test data
        sc.pp.normalize_total(adata_test, target_sum=1e4)
        sc.pp.log1p(adata_test)
        
        print("✓ Data normalized")

        # Start memory tracking
        tracemalloc.start()

        # Train CellTypist model
        print("Training CellTypist model...")
        print("Parameters: use_SGD=True, feature_selection=True")

        write_status(status_file, "TRAINING")

        # Train with specified parameters
        custom_model = celltypist.train(
            adata_train,
            labels='cell_type',
            use_SGD=True,
            feature_selection=True,
            n_jobs=4
        )

        print("✓ Model training completed")

        # Run predictions
        print("Running predictions...")

        write_status(status_file, "PREDICTING")

        predictions = celltypist.annotate(
            adata_test,
            model=custom_model,
            majority_voting=True
        )

        print("✓ Predictions completed")

        # Get peak memory usage
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_memory_mb = peak / (1024 * 1024)
        print(f"Peak memory usage: {peak_memory_mb:.2f} MB")
        print(f"PEAK_MEMORY_MB:{peak_memory_mb:.2f}")
        
        # Extract results
        print("Extracting results...")
        
        # Get AnnData with inserted labels and confidence scores
        adata_result = predictions.to_adata(insert_labels=True, insert_conf=True)
        
        # Extract predicted labels and confidence scores
        predicted_labels = adata_result.obs['predicted_labels'].values
        conf_scores = adata_result.obs['conf_score'].values
        cell_ids = adata_result.obs_names.values
        
        print(f"Extracted {len(predicted_labels)} predictions")
        print(f"Unique predicted cell types: {len(set(predicted_labels))}")
        print(f"Confidence score range: {np.min(conf_scores):.3f} - {np.max(conf_scores):.3f}")
        
        # Create results dataframe
        results_df = pd.DataFrame({
            'cell_id': cell_ids,
            'predicted_label': predicted_labels,
            'confidence_score': conf_scores
        })
        
        # Write results to CSV
        print(f"Writing results to: {output_file}")
        results_df.to_csv(output_file, index=False)
        
        print("✓ Results written successfully")
        
        # Write success status
        write_status(status_file, "SUCCESS")
        
        print("🎉 CellTypist processing completed successfully!")
        
    except Exception as e:
        error_msg = f"ERROR: {str(e)}"
        print(error_msg)
        write_status(status_file, error_msg)
        sys.exit(1)

def write_status(status_file, status_msg):
    """Write status message to status file"""
    try:
        with open(status_file, 'w') as f:
            f.write(status_msg + '\n')
    except Exception as e:
        print(f"Warning: Could not write status file: {e}")

if __name__ == "__main__":
    main()
