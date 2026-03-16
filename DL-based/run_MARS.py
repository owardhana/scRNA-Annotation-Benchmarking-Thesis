# run_MARS.py
#################################################
# MARS Function for Python Benchmarking Framework
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


def run_MARS_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    MARS Cell Type Annotation Function

    Purpose: Run MARS (Meta-learning Approach for Recognition in Single-cell data)
    Inputs:
      - adata_train: Training AnnData object (used for meta-learning)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses MARS's meta-learning approach to discover and annotate cell types
    Strengths: Discovers novel cell types, transfers knowledge across heterogeneous experiments
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
            "conda info --envs | grep MARS_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("MARS_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for MARS...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for MARS analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # MARS requires cell_type annotation in obs for labeled experiments
            # Ensure it's named 'cell_type' for MARS compatibility
            if 'cell_type' not in adata_train_subset.obs.columns:
                adata_train_subset.obs['cell_type'] = adata_train_subset.obs['Ground_Truth_Celltype']

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running MARS in conda environment...")

            # Create Python script for MARS execution
            mars_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    from mars import MARS

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['cell_type'].unique())}}")

    # Prepare MARS format
    # MARS expects:
    # - labeled_exp: List[AnnData] (labeled experiments)
    # - unlabeled_exp: AnnData (unlabeled target)
    # - pretrain_data: List[AnnData] or None

    print("Preparing MARS data format...")
    labeled_exp = [adata_train]  # Wrap train as single labeled experiment
    unlabeled_exp = adata_test    # Test as unlabeled target
    pretrain_data = [adata_train] # Use train as pretrain data

    # Auto-detect n_clusters from training data
    n_clusters = len(adata_train.obs['cell_type'].unique())
    print(f"Number of clusters (cell types): {{n_clusters}}")

    # MARS parameters
    # Note: params structure may need adjustment based on MARS version
    params = {{
        'use_cuda': False,  # Set to True if GPU available and desired
        'verbose': True
    }}

    print("Initializing MARS model...")
    # Initialize MARS
    try:
        mars = MARS(
            n_clusters=n_clusters,
            params=params,
            labeled_exp=labeled_exp,
            unlabeled_exp=unlabeled_exp,
            pretrain_data=pretrain_data
        )
    except TypeError as e:
        # MARS API may vary - try alternative initialization
        print(f"Standard initialization failed: {{e}}")
        print("Trying alternative MARS initialization...")
        mars = MARS(
            n_clusters=n_clusters,
            labeled_exp=labeled_exp,
            unlabeled_exp=unlabeled_exp,
            pretrain_data=pretrain_data
        )

    print("Training MARS model (this may take a while)...")
    # Train and predict (evaluation_mode=True restricts to known types)
    adata_result, landmarks, scores = mars.train(evaluation_mode=True)

    print("MARS training completed")

    # Extract predictions from adata_result
    # MARS may store predictions in various column names
    print("Extracting predictions...")
    print(f"Available columns in result: {{list(adata_result.obs.columns)}}")

    # Try to find prediction column
    prediction_columns = [
        'predicted_celltype', 'mars_prediction', 'celltype',
        'cell_type', 'predicted', 'pred_celltype'
    ]

    pred_col = None
    for col in prediction_columns:
        if col in adata_result.obs.columns:
            pred_col = col
            print(f"Found predictions in column: {{col}}")
            break

    if pred_col is None:
        # If no standard column found, try to use any column with 'pred' or 'type'
        for col in adata_result.obs.columns:
            if 'pred' in col.lower() or 'type' in col.lower():
                pred_col = col
                print(f"Using column {{col}} as predictions")
                break

    if pred_col is None:
        raise ValueError(f"Could not find prediction column. Available: {{list(adata_result.obs.columns)}}")

    predictions = list(adata_result.obs[pred_col])

    # Extract confidence scores if available
    # MARS may provide confidence in various ways
    confidence = None
    if scores is not None:
        # scores might be a matrix or dict
        if isinstance(scores, np.ndarray):
            # If scores is array, use max score per cell
            if scores.ndim == 2:
                confidence = np.max(scores, axis=1).tolist()
            else:
                confidence = scores.tolist()

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidence': confidence,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index)
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("MARS execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"MARS execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute MARS in conda environment
            result = subprocess.run(
                f'conda run -n MARS_env python -c "{mars_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=3600  # 60 minute timeout for meta-learning
            )

            if result.returncode != 0:
                warnings.warn(f"MARS execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing MARS results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("MARS results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidence_from_mars = results.get('confidence', None)
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Process confidence scores
            if confidence_from_mars is not None and len(confidence_from_mars) == len(cell_predictions):
                # Use MARS-provided confidence if available
                confidence_scores = []
                for conf in confidence_from_mars:
                    if pd.notna(conf) and isinstance(conf, (int, float)):
                        # Ensure confidence is between 0 and 1
                        conf_val = max(0.0, min(1.0, float(conf)))
                        confidence_scores.append(conf_val)
                    else:
                        confidence_scores.append(0.5)  # Default confidence
            else:
                # Use heuristic if MARS doesn't provide confidence
                confidence_scores = [0.8 if pred != "Unknown" else 0.0
                                   for pred in cell_predictions]

            # Adjust confidence for Unknown predictions
            for i, pred in enumerate(cell_predictions):
                if pred == "Unknown":
                    confidence_scores[i] = 0.0

            print(f"MARS completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            if confidence_scores:
                print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("MARS execution timed out after 60 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"MARS error: {str(e)}")
            return default_return()


# For backward compatibility
run_MARS = run_MARS_function
