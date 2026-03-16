# run_ItClust.py
#################################################
# ItClust Function for Python Benchmarking Framework
# Input: Train/test AnnData objects and markers from rank_genes_groups
# Output: Standardized results format for CV framework
#################################################

import subprocess
import tempfile
import os
import pandas as pd
import anndata as ad
import warnings
import pickle
import h5py
from typing import Dict, List, Any


def run_ItClust_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    ItClust Cell Type Annotation Function

    Purpose: Run ItClust algorithm using transfer learning with deep embedded clustering
    Inputs:
      - adata_train: Training AnnData object (used as source/reference data)
      - adata_test: Test AnnData object to predict (target data)
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses ItClust's transfer learning classifier with DEC (Deep Embedded Clustering)
               for transferring cell type annotations from source to target data
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
    conda_env_path = "/home/oliver/miniconda3/envs/ItClust_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"ItClust_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # ItClust requires "cell type" column for source data (not "Ground_Truth_Celltype")
            adata_train_subset.obs['cell type'] = adata_train_subset.obs['Ground_Truth_Celltype'].copy()

            # Save AnnData objects for inter-process communication
            # Strip to minimal data for old anndata (0.8.0) compatibility
            adata_train_subset.obs = adata_train_subset.obs[['Ground_Truth_Celltype', 'cell type']].copy()
            adata_test_subset.obs = adata_test_subset.obs[['Ground_Truth_Celltype']].copy()
            for adata in [adata_train_subset, adata_test_subset]:
                for col in adata.obs.columns:
                    adata.obs[col] = adata.obs[col].astype(str)
                adata.var = adata.var[[]]
                adata.uns.clear()
                adata.layers.clear()
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)
            # Remove /layers group from h5ad — old anndata can't parse encoding-type 'dict'
            for h5_path in [train_h5ad_path, test_h5ad_path]:
                with h5py.File(h5_path, 'a') as f:
                    if 'layers' in f:
                        del f['layers']

            # Create Python script for ItClust execution
            itclust_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory

import sys
import os
import pickle
import warnings
warnings.filterwarnings("ignore")

# Set seeds for reproducibility
from numpy.random import seed
seed(20180806)
import numpy as np
np.random.seed(10)

try:
    import tensorflow as tf
    # For TensorFlow 1.x compatibility
    if hasattr(tf, 'set_random_seed'):
        tf.set_random_seed(20180806)
    else:
        tf.random.set_seed(20180806)

    import ItClust as ic
    import scanpy as sc
    import anndata as ad
    import pandas as pd

    # Load data
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # ItClust requires obs["celltype"] as the label column
    adata_train.obs['celltype'] = adata_train.obs['Ground_Truth_Celltype']

    # ============================================
    # Fit ItClust model and get predictions
    # ============================================
    clf = ic.transfer_learning_clf()
    clf.fit(adata_train, adata_test)
    pred, prob, cell_type_pred = clf.predict()

    # cell_type_pred is dict: {{cluster_id_str: [celltype, confidence]}}
    # pred['cluster'] holds each test cell's assigned cluster id (as string)
    cluster_ids = pred['cluster'].tolist()
    predictions = [cell_type_pred.get(c, ['Unknown', 0.0])[0] for c in cluster_ids]
    confidence_scores = [cell_type_pred.get(c, ['Unknown', 0.0])[1] for c in cluster_ids]
    cell_ids = pred['cell_id'].tolist()
    true_labels = list(adata_test.obs['Ground_Truth_Celltype'])

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidence_scores': confidence_scores,
        'true_labels': true_labels,
        'cell_ids': cell_ids
    }}

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)  # Convert bytes to MB

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print(f"Predictions: {{len(predictions)}} cells, unique: {{set(predictions)}}")

except Exception as e:
    print(f"ItClust execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute ItClust in conda environment
            # Set environment - enable GPU
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend
            # Note: ItClust uses TensorFlow, GPU will be used if available

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "ItClust_env",
                    "python",
                    "-c", itclust_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"ItClust execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("ItClust results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidence_scores = results['confidence_scores']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            peak_memory_mb = results.get('peak_memory_mb', None)

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Process confidence scores
            confidence_scores = [
                0.0 if pred == "Unknown" else (max(0.0, min(1.0, float(conf))) if pd.notna(conf) and isinstance(conf, (int, float)) else 0.5)
                for pred, conf in zip(cell_predictions, confidence_scores)
            ]

            print(f"ItClust completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb
            }

        except subprocess.TimeoutExpired:
            warnings.warn("ItClust execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"ItClust error: {str(e)}")
            return default_return()


# For backward compatibility
run_ItClust = run_ItClust_function
