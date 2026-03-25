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
import h5py
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

    # Check conda environment availability (direct directory check)
    conda_env_path = "/home/oliver/miniconda3/envs/scBalance_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scBalance_env conda environment not found at {conda_env_path}. Please create it first.")
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
            # Strip to minimal data for old anndata (0.10.9) compatibility
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

            print("Running scBalance in conda environment...")

            # Create Python script for scBalance execution
            scBalance_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle
import time
import torch
from memory_profiler import memory_usage

try:
    import scBalance as sb
    import scipy.sparse as sp

    # Load data (NOT measured)
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # Ensure float32 dtype (NOT measured)
    for adata in [adata_train, adata_test]:
        if sp.issparse(adata.X):
            adata.X = adata.X.astype(np.float32)
        else:
            adata.X = np.array(adata.X, dtype=np.float32)

    # Build DataFrames for scBalance (NOT measured)
    train_X = np.array(adata_train.X.todense()) if hasattr(adata_train.X, 'todense') else np.array(adata_train.X)
    test_X  = np.array(adata_test.X.todense())  if hasattr(adata_test.X,  'todense') else np.array(adata_test.X)
    reference = pd.DataFrame(train_X, index=adata_train.obs.index, columns=adata_train.var.index)
    test      = pd.DataFrame(test_X,  index=adata_test.obs.index,  columns=adata_test.var.index)
    celltype_values = adata_train.obs['Ground_Truth_Celltype']
    celltype_array  = celltype_values.cat.categories[celltype_values.cat.codes].values if hasattr(celltype_values, 'cat') else celltype_values.values
    label = pd.DataFrame({{'cell_type': celltype_array}}, index=adata_train.obs.index)

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        processing_unit = 'gpu' if torch.cuda.is_available() else 'cpu'
        pred_result = sb.scBalance(
            test=test, reference=reference, label=label,
            processing_unit=processing_unit, save_model=False, weighted_sampling=True)

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return list(pred_result), list(adata_test.obs['Ground_Truth_Celltype']), list(adata_test.obs.index), _peak_vram_mb

    mem_list, (pred_result, true_labels, cell_ids, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': pred_result,
        'true_labels': true_labels,
        'cell_ids': cell_ids,
        'peak_system_memory_mb': peak_system_memory_mb,
        'peak_vram_mb': peak_vram_mb,
        'method_walltime': method_walltime,
    }}

    print(f"Peak system memory: {{peak_system_memory_mb:.2f}} MB, VRAM: {{peak_vram_mb:.2f}} MB, walltime: {{method_walltime:.2f}}s")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print(f"Predictions: {{len(results['predictions'])}} cells, unique: {{set(results['predictions'])}}")

except Exception as e:
    print(f"scBalance execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scBalance in conda environment
            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scBalance_env",
                    "python",
                    "-c", scBalance_script
                ],
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scBalance execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

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
                'cell_ids': cell_ids,
                'peak_memory_mb':  results.get('peak_system_memory_mb', None),
                'peak_vram_mb':    results.get('peak_vram_mb', None),
                'method_walltime': results.get('method_walltime', None),
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scBalance execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scBalance error: {str(e)}")
            return default_return()


# For backward compatibility
run_scBalance = run_scBalance_function
