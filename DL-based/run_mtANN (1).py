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
import warnings
import pickle
import h5py
from typing import Dict, List, Any


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

    # Check conda environment availability (direct directory check)
    conda_env_path = "/home/oliver/miniconda3/envs/mtANN_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"mtANN_env conda environment not found at {conda_env_path}. Please create it first.")
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

            # Save AnnData objects for inter-process communication
            # Strip to minimal data for old anndata (0.9.2) compatibility
            for adata in [adata_train_subset, adata_test_subset]:
                adata.obs = adata.obs[['Ground_Truth_Celltype']].copy()
                adata.obs['Ground_Truth_Celltype'] = adata.obs['Ground_Truth_Celltype'].astype(str)
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

            # Create Python script for mtANN execution
            mtann_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory (before imports, after data conversion)

import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import scanpy as sc
    from mtANN.model import mtANN

    # Load data
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # Check for CUDA availability
    try:
        import torch
        cuda_available = torch.cuda.is_available()
    except ImportError:
        cuda_available = False

    # Prepare expression matrices for mtANN

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

    # Clean data to remove NaN, Inf, and extreme values
    train_expr = np.nan_to_num(train_expr, nan=0.0, posinf=0.0, neginf=0.0)
    train_expr = np.clip(train_expr, -1e10, 1e10)

    test_expr = np.nan_to_num(test_expr, nan=0.0, posinf=0.0, neginf=0.0)
    test_expr = np.clip(test_expr, -1e10, 1e10)

    # mtANN expects pandas DataFrames with gene names as columns
    import pandas as pd
    gene_names = adata_train.var.index.tolist()

    train_df = pd.DataFrame(train_expr,
                           index=adata_train.obs.index,
                           columns=gene_names)
    test_df = pd.DataFrame(test_expr,
                          index=adata_test.obs.index,
                          columns=gene_names)

    # Prepare mtANN input format
    # mtANN's gene_select="default" requires an R script at a hardcoded path that doesn't exist.
    # gene_select="manual" with 1 reference creates only 1 classifier instead of the 8 mtANN needs.
    # Workaround: create 8 diverse gene subsets and pass as separate "references" for ensemble diversity.
    labels = adata_train.obs["Ground_Truth_Celltype"].values
    n_genes = len(gene_names)
    subset_size = min(n_genes, 2000)

    gene_subsets = []
    # 1. Top genes by variance
    gene_subsets.append(train_df.var().nlargest(subset_size).index.tolist())
    # 2. Top genes by mean expression
    gene_subsets.append(train_df.mean().nlargest(subset_size).index.tolist())
    # 3. Top genes by coefficient of variation
    cv = train_df.std() / (train_df.mean() + 1e-8)
    gene_subsets.append(cv.nlargest(subset_size).index.tolist())
    # 4. Top genes by max expression
    gene_subsets.append(train_df.max().nlargest(subset_size).index.tolist())
    # 5-8. Random subsets with different seeds
    for seed_val in [42, 123, 456, 789]:
        rng = np.random.RandomState(seed_val)
        gene_subsets.append(rng.choice(gene_names, size=subset_size, replace=False).tolist())

    expression_s = [train_df[subset] for subset in gene_subsets]
    label_s = [labels] * len(gene_subsets)
    expression_t = test_df  # mtANN intersects with each subset's columns

    # Run mtANN with 8-classifier ensemble
    mid_annotation, final_annotation, m, threshold = mtANN(
        expression_s=expression_s,
        label_s=label_s,
        expression_t=expression_t,
        threshold=0.0,  # Fixed threshold to bypass buggy threshold_selection
        gene_select="manual",  # Each reference has different genes, creating 8 diverse classifiers
        CUDA=cuda_available
    )

    # Process results (use final_annotation if available, otherwise mid_annotation)
    predictions = final_annotation if final_annotation is not None else mid_annotation
    if predictions is None:
        raise ValueError("No annotations returned from mtANN")

    # Convert to list and flatten if needed
    if isinstance(predictions, np.ndarray):
        predictions = predictions.flatten().tolist()

    # Ensure predictions are strings (not arrays or lists)
    predictions = [str(p) if not isinstance(p, str) else p for p in predictions]

    # Handle nested lists if they exist
    if predictions and isinstance(predictions[0], (list, np.ndarray)):
        predictions = [str(p[0]) if len(p) > 0 else "Unknown" for p in predictions]

    # Force final conversion to clean strings (strip brackets and quotes from string representation)
    predictions = [str(p).strip('[]').strip("'").strip('"') for p in predictions]

    # Prepare results
    results = {{
        'predictions': predictions,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index),
        'threshold': threshold
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
    print(f"mtANN execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute mtANN in conda environment
            # Set environment to avoid Jupyter backend conflicts
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "mtANN_env",
                    "python",
                    "-c", mtann_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout
            )

            if result.returncode != 0:
                warnings.warn(f"mtANN execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

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
            peak_memory_mb = results.get('peak_memory_mb', None)  # Extract peak memory

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores (mtANN doesn't provide explicit confidence)
            confidence_scores = [0.8 if pred != "Unknown" else 0.0 for pred in cell_predictions]

            print(f"mtANN completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb  # Add peak memory to return dict
            }

        except subprocess.TimeoutExpired:
            warnings.warn("mtANN execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"mtANN error: {str(e)}")
            return default_return()


# For backward compatibility
run_mtANN = run_mtANN_function