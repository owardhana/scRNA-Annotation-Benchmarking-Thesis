# run_CIForm.py
#################################################
# CIForm Function for Python Benchmarking Framework
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


def run_CIForm_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    CIForm Cell Type Annotation Function

    Purpose: Run CIForm algorithm using consensus integration with AnnData-compatible helper
    Inputs:
      - adata_train: Training AnnData object (used as reference dataset)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses CIForm's consensus integration approach with AnnData objects
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
            "conda info --envs | grep CIForm_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("CIForm_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for CIForm...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for CIForm analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running CIForm in conda environment...")

            # Create Python script for CIForm execution using the helper
            ciform_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    # Import the CIForm helper with AnnData-compatible functions
    from run_CIForm_helper import ciForm

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Prepare training labels
    train_labels = adata_train.obs['Ground_Truth_Celltype'].tolist()

    print("Running CIForm consensus integration...")

    # Run CIForm with AnnData-compatible function
    pred_result = ciForm(
        s=1024,  # Length of sub-vector
        Train_adata=adata_train,
        train_labels=train_labels,
        Test_adata=adata_test,
        n_epochs=20
    )

    print("CIForm execution completed")

    # Process results
    if isinstance(pred_result, (list, np.ndarray)):
        predictions = list(pred_result)
    elif isinstance(pred_result, dict):
        # If result is a dictionary, try to extract predictions
        if 'predictions' in pred_result:
            predictions = list(pred_result['predictions'])
        elif 'labels' in pred_result:
            predictions = list(pred_result['labels'])
        else:
            # Use first value if structure is unknown
            predictions = list(list(pred_result.values())[0])
    else:
        # Handle single prediction case
        predictions = [str(pred_result)]

    # Ensure prediction length matches test data
    if len(predictions) != len(adata_test.obs):
        # If mismatch, pad or truncate predictions
        if len(predictions) < len(adata_test.obs):
            predictions.extend(['Unknown'] * (len(adata_test.obs) - len(predictions)))
        else:
            predictions = predictions[:len(adata_test.obs)]

    # Prepare results
    results = {{
        'predictions': predictions,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index)
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("CIForm execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"CIForm execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute CIForm in conda environment
            result = subprocess.run(
                f'conda run -n CIForm_env python -c "{ciform_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minute timeout
            )

            if result.returncode != 0:
                warnings.warn(f"CIForm execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing CIForm results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("CIForm results file not found")
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

            # Create confidence scores (CIForm doesn't provide explicit confidence)
            # Use heuristic based on consensus prediction quality
            confidence_scores = []
            training_cell_types = set(adata_train.obs['Ground_Truth_Celltype'].unique())

            for pred in cell_predictions:
                if pred == "Unknown":
                    confidence_scores.append(0.0)
                elif pred in training_cell_types:
                    confidence_scores.append(0.7)  # Moderate confidence for consensus predictions
                else:
                    confidence_scores.append(0.4)  # Lower confidence for novel types

            print(f"CIForm completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("CIForm execution timed out after 30 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"CIForm error: {str(e)}")
            return default_return()


# For backward compatibility
run_CIForm = run_CIForm_function