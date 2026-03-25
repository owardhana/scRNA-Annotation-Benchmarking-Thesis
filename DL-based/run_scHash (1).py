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
import anndata as ad
import warnings
import pickle
from typing import Dict, List, Any


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

    # Check conda environment availability (direct directory check)
    conda_env_path = "/home/oliver/miniconda3/envs/scHash_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scHash_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
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

            # Subset to common genes with consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Add dummy batch column (scHash requires batch_key even if no batches exist)
            adata_train_subset.obs['batch'] = 'batch1'
            adata_test_subset.obs['batch'] = 'batch1'

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            # Create Python script for scHash execution
            schash_script = f'''
import sys
import pickle
import time
import torch
from memory_profiler import memory_usage

try:
    import scHash
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    # Load data (NOT measured)
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # Ensure float32 dtype (scHash's StandardScaler outputs float64, but PyTorch Linear layers expect float32)
    for adata in [adata_train, adata_test]:
        if sp.issparse(adata.X):
            adata.X = adata.X.astype(np.float32)
        else:
            adata.X = np.array(adata.X, dtype=np.float32)

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        datamodule = scHash.setup_training_data(
            train_data=adata_train,
            cell_type_key='Ground_Truth_Celltype',
            batch_key='batch',
            log_norm=False,
            hvg=False,
            normalize=True,
            batch_size=128
        )
        model = scHash.scHashModel(datamodule)
        trainer, best_model_path, training_time = scHash.training(
            model=model, datamodule=datamodule,
            checkpointPath='{checkpoint_path}', max_epochs=50)
        datamodule.setup_test_data(adata_test)
        pred_labels, hash_codes = scHash.testing(trainer, model, best_model_path)

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return pred_labels, hash_codes, list(adata_test.obs['Ground_Truth_Celltype']), list(adata_test.obs.index), training_time, _peak_vram_mb

    mem_list, (pred_labels, hash_codes, true_labels, cell_ids, training_time, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': pred_labels,
        'hash_codes': hash_codes,
        'true_labels': true_labels,
        'cell_ids': cell_ids,
        'training_time': training_time,
        'peak_system_memory_mb': peak_system_memory_mb,
        'peak_vram_mb': peak_vram_mb,
        'method_walltime': method_walltime,
    }}

    print(f"Peak system memory: {{peak_system_memory_mb:.2f}} MB, VRAM: {{peak_vram_mb:.2f}} MB, walltime: {{method_walltime:.2f}}s")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print(f"Predictions: {{len(pred_labels)}} cells, unique: {{set(pred_labels)}}")

except Exception as e:
    print(f"scHash execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scHash in conda environment
            # Set environment to avoid Jupyter backend conflicts
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scHash_env",
                    "python",
                    "-c", schash_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scHash execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

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
            confidence_scores = [0.8 if pred != "Unknown" else 0.0 for pred in cell_predictions]

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
                'cell_ids': cell_ids,
                'peak_memory_mb':  results.get('peak_system_memory_mb', None),
                'peak_vram_mb':    results.get('peak_vram_mb', None),
                'method_walltime': results.get('method_walltime', None),
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scHash execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scHash error: {str(e)}")
            return default_return()


# For backward compatibility
run_scHash = run_scHash_function