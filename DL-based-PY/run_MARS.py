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

    # Check conda environment availability (direct directory check)
    conda_env_path = "/home/oliver/miniconda3/envs/MARS_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"MARS_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for MARS...")

            # Prepare file paths (use npy + json for maximum compatibility)
            train_X_path = os.path.join(temp_dir, "train_X.npy")
            train_meta_path = os.path.join(temp_dir, "train_meta.json")
            test_X_path = os.path.join(temp_dir, "test_X.npy")
            test_meta_path = os.path.join(temp_dir, "test_meta.json")
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

            # Save data using npy + json for maximum compatibility
            # Avoids h5ad, pickle, and numpy version issues
            import json

            # Get expression matrices as plain numpy arrays
            if hasattr(adata_train_subset.X, 'todense'):
                train_X = np.array(adata_train_subset.X.todense(), dtype=np.float32)
            else:
                train_X = np.array(adata_train_subset.X, dtype=np.float32)

            if hasattr(adata_test_subset.X, 'todense'):
                test_X = np.array(adata_test_subset.X.todense(), dtype=np.float32)
            else:
                test_X = np.array(adata_test_subset.X, dtype=np.float32)

            # Save train data - expression matrix as .npy, metadata as JSON
            np.save(train_X_path, train_X)
            with open(train_meta_path, 'w') as f:
                json.dump({
                    'cell_types': adata_train_subset.obs['cell_type'].astype(str).tolist(),
                    'cell_ids': adata_train_subset.obs.index.astype(str).tolist(),
                    'gene_ids': adata_train_subset.var.index.astype(str).tolist()
                }, f)

            # Save test data
            np.save(test_X_path, test_X)
            with open(test_meta_path, 'w') as f:
                json.dump({
                    'cell_ids': adata_test_subset.obs.index.astype(str).tolist(),
                    'gene_ids': adata_test_subset.var.index.astype(str).tolist()
                }, f)

            print("Running MARS in conda environment...")

            # Create Python script for MARS execution
            mars_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory (before imports, after data conversion)

import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

# Add MARS package to Python path
mars_package_path = '/home/oliver/Thesis/benchmarking/DL-based/MARS_package'
if mars_package_path not in sys.path:
    sys.path.insert(0, mars_package_path)

try:
    # Import MARS from cloned repository
    from model.mars import MARS
    from model.experiment_dataset import ExperimentDataset

    # Load data from npy + json files
    print("Loading training and test data...")
    import json

    train_X = np.load('{train_X_path}')
    with open('{train_meta_path}', 'r') as f:
        train_meta = json.load(f)
    train_cell_types = train_meta['cell_types']
    train_cell_ids = train_meta['cell_ids']
    train_gene_ids = train_meta['gene_ids']

    test_X = np.load('{test_X_path}')
    with open('{test_meta_path}', 'r') as f:
        test_meta = json.load(f)
    test_cell_ids = test_meta['cell_ids']
    test_gene_ids = test_meta['gene_ids']

    print(f"Training data: {{train_X.shape[0]}} cells x {{train_X.shape[1]}} genes")
    print(f"Test data: {{test_X.shape[0]}} cells x {{test_X.shape[1]}} genes")
    print(f"Training cell types: {{list(np.unique(train_cell_types))}}")

    # Prepare MARS format
    # MARS expects:
    # - labeled_data: List[ExperimentDataset] (labeled experiments)
    # - unlabeled_data: ExperimentDataset (unlabeled target)
    # - pretrain_data: List[ExperimentDataset] or None

    print("Preparing MARS data format...")

    # Convert to ExperimentDataset format
    # ExperimentDataset(x, cells, genes, metadata, y=[])
    # x: rows=cells, cols=genes
    # y: numeric labels for labeled data

    # Convert cell types to numeric labels
    cell_types = np.unique(train_cell_types)
    cell_type_to_idx = {{ct: idx for idx, ct in enumerate(cell_types)}}
    idx_to_cell_type = {{idx: ct for ct, idx in cell_type_to_idx.items()}}
    y_train = np.array([cell_type_to_idx[ct] for ct in train_cell_types])

    # Create ExperimentDataset objects
    train_dataset = ExperimentDataset(
        x=train_X,
        cells=list(train_cell_ids),
        genes=list(train_gene_ids),
        metadata='train',
        y=y_train
    )

    test_dataset = ExperimentDataset(
        x=test_X,
        cells=list(test_cell_ids),
        genes=list(test_gene_ids),
        metadata='test',
        y=[]  # No labels for test data
    )

    labeled_data = [train_dataset]  # Wrap train as single labeled experiment
    unlabeled_data = test_dataset    # Test as unlabeled target
    pretrain_data = [train_dataset] # Use train as pretrain data

    # Auto-detect n_clusters from training data
    n_clusters = len(cell_types)
    print(f"Number of clusters (cell types): {{n_clusters}}")

    # MARS parameters - must be object with attributes (not dict)
    import argparse
    import torch
    params = argparse.Namespace()
    params.pretrain_batch = None  # No batching for pretraining
    params.device = torch.device('cpu')  # Use CPU (RTX 5060 sm_120 incompatible with PyTorch 1.2.0)
    params.epochs = 30  # Training epochs
    params.epochs_pretrain = 25  # Pretraining epochs
    params.pretrain = True  # Enable pretraining
    params.model_file = '/tmp/mars_model.pt'  # Temporary model file
    params.learning_rate = 0.001  # Learning rate
    params.lr_scheduler_gamma = 0.5  # LR scheduler gamma
    params.lr_scheduler_step = 20  # LR scheduler step

    print("Initializing MARS model...")
    # Initialize MARS
    # MARS signature: __init__(self, n_clusters, params, labeled_data, unlabeled_data, pretrain_data=None, ...)
    mars = MARS(
        n_clusters=n_clusters,
        params=params,
        labeled_data=labeled_data,
        unlabeled_data=unlabeled_data,
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

    predictions_raw = list(adata_result.obs[pred_col])

    # Convert numeric predictions to cell type names if needed
    predictions = []
    for pred in predictions_raw:
        if isinstance(pred, (int, np.integer)):
            # Numeric prediction - map to cell type
            predictions.append(idx_to_cell_type.get(pred, "Unknown"))
        else:
            # Already a string
            predictions.append(str(pred))

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

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)  # Convert bytes to MB

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

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
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "MARS_env",
                    "python",
                    "-c", mars_script
                ],
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout for meta-learning
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
            peak_memory_mb = results.get('peak_memory_mb', None)  # Extract peak memory

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
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb  # Add peak memory to return dict
            }

        except subprocess.TimeoutExpired:
            warnings.warn("MARS execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"MARS error: {str(e)}")
            return default_return()


# For backward compatibility
run_MARS = run_MARS_function
