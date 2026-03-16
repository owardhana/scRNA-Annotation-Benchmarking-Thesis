# run_scHash.py
#################################################
# scHash Function for Python Benchmarking Framework
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


def run_scHash_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scHash Cell Type Annotation Function

    Purpose: Run scHash algorithm using hash-based deep learning with lightning framework
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scHash's hash-based deep learning classifier with conda environment activation
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
            "conda info --envs | grep scHash_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("scHash_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scHash...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            checkpoint_path = os.path.join(temp_dir, "checkpoint")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Create checkpoint directory
            os.makedirs(checkpoint_path, exist_ok=True)

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for scHash analysis")

            # Subset both datasets to common genes
            adata_train_subset = adata_train[:, common_genes].copy()
            adata_test_subset = adata_test[:, common_genes].copy()

            # Ensure gene order is consistent
            adata_train_subset = adata_train_subset[:, sorted(common_genes)]
            adata_test_subset = adata_test_subset[:, sorted(common_genes)]

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running scHash in conda environment...")

            # Create Python script for scHash execution
            schash_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle
from statistics import median

try:
    import scHash
    from sklearn.metrics import f1_score, precision_score, recall_score

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Set up the training datamodule
    print("Setting up scHash training datamodule...")
    datamodule = scHash.setup_training_data(
        train_data=adata_train,
        cell_type_key='Ground_Truth_Celltype'
    )

    # Initialize scHash model
    print("Initializing scHash model...")
    model = scHash.scHashModel(datamodule)

    # Train the model with reduced epochs for cross-validation efficiency
    print("Training scHash model...")
    trainer, best_model_path, training_time = scHash.training(
        model=model,
        datamodule=datamodule,
        checkpointPath='{checkpoint_path}',
        max_epochs=50  # Reduced for CV efficiency
    )

    # Add the test data to datamodule
    print("Setting up test data...")
    datamodule.setup_test_data(adata_test)

    # Test the model
    print("Making predictions...")
    pred_labels, hash_codes = scHash.testing(trainer, model, best_model_path)

    # Prepare results
    results = {{
        'predictions': pred_labels,
        'hash_codes': hash_codes,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index),
        'training_time': training_time
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("scHash execution completed successfully")
    print(f"Predictions: {{len(pred_labels)}} cells")
    print(f"Unique predictions: {{set(pred_labels)}}")

except Exception as e:
    print(f"scHash execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scHash in conda environment
            result = subprocess.run(
                f'conda run -n scHash_env python -c "{schash_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=2400  # 40 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scHash execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing scHash results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scHash results file not found")
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

            # Create confidence scores (scHash doesn't provide explicit confidence)
            # Use heuristic based on training cell types
            training_cell_types = set(adata_train.obs['Ground_Truth_Celltype'].unique())
            confidence_scores = []

            for pred in cell_predictions:
                if pred == "Unknown":
                    confidence_scores.append(0.0)
                elif pred in training_cell_types:
                    confidence_scores.append(0.8)  # High confidence for known types
                else:
                    confidence_scores.append(0.3)  # Lower confidence for novel types

            print(f"scHash completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Training time: {results.get('training_time', 'N/A')} seconds")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scHash execution timed out after 40 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"scHash error: {str(e)}")
            return default_return()


# For backward compatibility
run_scHash = run_scHash_function