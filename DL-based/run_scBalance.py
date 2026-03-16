# run_scBalance.py
#################################################
# scBalance Function for Python Benchmarking Framework
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


def run_scBalance_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scBalance Cell Type Annotation Function

    Purpose: Run scBalance algorithm using dropout neural network for cell type classification
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scBalance's dropout neural network with weighted sampling for rare cell types
    Strengths: Excellent performance on rare cell types and batch effect robustness
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
            "conda info --envs | grep scBalance_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("scBalance_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scBalance...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for scBalance analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running scBalance in conda environment...")

            # Create Python script for scBalance execution
            scBalance_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import scBalance as sb
    import torch

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Convert AnnData to pandas DataFrames (scBalance format)
    print("Converting AnnData to pandas format...")

    # scBalance expects:
    # - reference: genes (rows) x cells (columns)
    # - test: genes (rows) x cells (columns)
    # - label: pandas Series/DataFrame with cell type labels

    # Get expression matrices
    # Handle both sparse and dense matrices
    if hasattr(adata_train.X, 'todense'):
        train_X = np.array(adata_train.X.todense())
    else:
        train_X = np.array(adata_train.X)

    if hasattr(adata_test.X, 'todense'):
        test_X = np.array(adata_test.X.todense())
    else:
        test_X = np.array(adata_test.X)

    # Create DataFrames: genes x cells format
    reference = pd.DataFrame(
        train_X.T,  # Transpose: cells x genes → genes x cells
        index=adata_train.var.index,
        columns=adata_train.obs.index
    )

    test = pd.DataFrame(
        test_X.T,  # Transpose: cells x genes → genes x cells
        index=adata_test.var.index,
        columns=adata_test.obs.index
    )

    # Create label DataFrame
    label = pd.DataFrame({{
        'cell_type': adata_train.obs['Ground_Truth_Celltype'].values
    }}, index=adata_train.obs.index)

    print(f"Reference shape: {{reference.shape}}")
    print(f"Test shape: {{test.shape}}")
    print(f"Label shape: {{label.shape}}")

    # Auto-detect processing unit (GPU if available, else CPU)
    processing_unit = 'gpu' if torch.cuda.is_available() else 'cpu'
    print(f"Using processing unit: {{processing_unit}}")
    if processing_unit == 'gpu':
        print(f"GPU: {{torch.cuda.get_device_name(0)}}")

    # Run scBalance
    print("Training scBalance model...")
    pred_result = sb.scBalance(
        test=test,
        reference=reference,
        label=label,
        processing_unit=processing_unit,
        save_model=False,  # Don't save model for benchmarking
        weighted_sampling=True  # Handle class imbalance (good for rare types)
    )

    print("scBalance prediction completed")

    # Prepare results
    # pred_result should be array-like with predictions
    results = {{
        'predictions': list(pred_result),
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index)
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("scBalance execution completed successfully")
    print(f"Predictions: {{len(results['predictions'])}} cells")
    print(f"Unique predictions: {{set(results['predictions'])}}")

except Exception as e:
    print(f"scBalance execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scBalance in conda environment
            result = subprocess.run(
                f'conda run -n scBalance_env python -c "{scBalance_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=3600  # 60 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scBalance execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing scBalance results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scBalance results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # scBalance doesn't provide confidence scores
            # Use heuristic: higher confidence for valid predictions
            confidence_scores = []
            for pred in cell_predictions:
                if pred == "Unknown":
                    confidence_scores.append(0.0)
                else:
                    # Default confidence of 0.8 for scBalance predictions
                    # (scBalance is known for high accuracy, especially on rare types)
                    confidence_scores.append(0.8)

            print(f"scBalance completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scBalance execution timed out after 60 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"scBalance error: {str(e)}")
            return default_return()


# For backward compatibility
run_scBalance = run_scBalance_function
