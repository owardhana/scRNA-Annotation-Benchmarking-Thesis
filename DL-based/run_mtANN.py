# run_mtANN.py
#################################################
# mtANN Function for Python Benchmarking Framework
# Input: Train/test AnnData objects and markers from rank_genes_groups
# Output: Standardized results format for CV framework
#################################################

import subprocess
import tempfile
import os
import pandas as pd
import numpy as np
import anndata as ad
import time
import warnings
import pickle
from typing import Dict, List, Any, Optional


def run_mtANN_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    mtANN Cell Type Annotation Function

    Purpose: Run mtANN algorithm using multi-task artificial neural network
    Inputs:
      - adata_train: Training AnnData object (used as reference for training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses mtANN's multi-task learning approach with expression matrices
    """

    # Default return function for error handling
    def default_return():
        n_test_cells = len(adata_test.obs)
        return {
            'predictions': ['Unknown'] * n_test_cells,
            'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
            'confidence_scores': [0.0] * n_test_cells,
            'cell_ids': list(adata_test.obs.index)
        }

    # Validate input data
    if 'Ground_Truth_Celltype' not in adata_train.obs.columns:
        warnings.warn("Ground_Truth_Celltype not found in training data")
        return default_return()

    if 'Ground_Truth_Celltype' not in adata_test.obs.columns:
        warnings.warn("Ground_Truth_Celltype not found in test data")
        return default_return()

    # Check conda environment availability
    try:
        env_check = subprocess.run(
            "conda info --envs | grep mtANN_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("mtANN_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for mtANN...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for mtANN analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running mtANN in conda environment...")

            # Create Python script for mtANN execution
            mtann_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import scanpy as sc
    from mtANN import mtANN

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Check for CUDA availability
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        print(f"CUDA available: {{cuda_available}}")
    except ImportError:
        cuda_available = False
        print("PyTorch not available, CUDA check skipped")

    # Prepare expression matrices for mtANN
    print("Preparing expression matrices...")

    # Convert to dense matrices if sparse
    if hasattr(adata_train.X, 'todense'):
        train_expr = adata_train.X.todense()
    else:
        train_expr = adata_train.X

    if hasattr(adata_test.X, 'todense'):
        test_expr = adata_test.X.todense()
    else:
        test_expr = adata_test.X

    # Convert to numpy arrays
    train_expr = np.array(train_expr)
    test_expr = np.array(test_expr)

    # Prepare mtANN input format
    expression_s = [train_expr]  # List of reference expression matrices
    label_s = [adata_train.obs["Ground_Truth_Celltype"].tolist()]  # List of label arrays
    expression_t = test_expr  # Target expression matrix

    print("Running mtANN annotation...")

    # Run mtANN with appropriate parameters
    mid_annotation, final_annotation, m, threshold = mtANN(
        expression_s=expression_s,
        label_s=label_s,
        expression_t=expression_t,
        threshold="default",
        gene_select="default",
        CUDA=cuda_available
    )

    print("mtANN execution completed")

    # Process results - use final_annotation as the primary result
    if final_annotation is not None:
        predictions = final_annotation
        annotation_type = "final"
    elif mid_annotation is not None:
        predictions = mid_annotation
        annotation_type = "mid"
    else:
        raise ValueError("No annotations returned from mtANN")

    print(f"Using {{annotation_type}} annotation results")

    # Ensure predictions is a list and has correct format
    if isinstance(predictions, np.ndarray):
        predictions = predictions.tolist()
    elif not isinstance(predictions, list):
        predictions = [predictions]

    # Handle case where predictions might be indices instead of labels
    if all(isinstance(p, (int, np.integer)) for p in predictions):
        # Convert indices to labels using training labels
        unique_labels = list(set(label_s[0]))
        predictions = [unique_labels[min(p, len(unique_labels)-1)] if p < len(unique_labels) else "Unknown"
                      for p in predictions]

    # Prepare results
    results = {{
        'predictions': predictions,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index),
        'threshold': threshold,
        'annotation_type': annotation_type
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("mtANN execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")
    print(f"Threshold used: {{threshold}}")

except Exception as e:
    print(f"mtANN execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute mtANN in conda environment
            result = subprocess.run(
                f'conda run -n mtANN_env python -c "{mtann_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minute timeout
            )

            if result.returncode != 0:
                warnings.warn(f"mtANN execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing mtANN results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("mtANN results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            threshold = results.get('threshold', 'default')
            annotation_type = results.get('annotation_type', 'unknown')

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores (mtANN doesn't provide explicit confidence)
            # Use heuristic based on annotation type and threshold
            confidence_scores = []
            training_cell_types = set(adata_train.obs['Ground_Truth_Celltype'].unique())

            base_confidence = 0.8 if annotation_type == "final" else 0.6

            for pred in cell_predictions:
                if pred == "Unknown":
                    confidence_scores.append(0.0)
                elif pred in training_cell_types:
                    confidence_scores.append(base_confidence)
                else:
                    # Novel cell type
                    confidence_scores.append(base_confidence * 0.7)

            print(f"mtANN completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Annotation type: {annotation_type}")
            print(f"  - Threshold: {threshold}")
            print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("mtANN execution timed out after 30 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"mtANN error: {str(e)}")
            return default_return()


# For backward compatibility
run_mtANN = run_mtANN_function