# run_scnym.py
#################################################
# scnym Function for Python Benchmarking Framework
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


def run_scnym_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scnym Cell Type Annotation Function

    Purpose: Run scnym algorithm using adversarial domain adaptation for novel cell type detection
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scnym's adversarial domain adaptation with novel cell type detection
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
    conda_env_path = "/home/oliver/miniconda3/envs/scnym_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scnym_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            combined_h5ad_path = os.path.join(temp_dir, "adata_combined.h5ad")
            temp_folder = os.path.join(temp_dir, "temp_folder")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Create temporary folder for scnym
            os.makedirs(temp_folder, exist_ok=True)

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Re-normalize to CPM (1e6) for scnym requirement
            # scnym requires log(CPM+1) specifically, not log(1e4+1)

            import numpy as np
            import scipy.sparse as sp

            for adata_subset in [adata_train_subset, adata_test_subset]:
                # Check if data is log-normalized (has log1p flag or max < 20)
                if adata_subset.X.max() < 20:  # Likely log-normalized
                    # Reverse log1p: exp(x) - 1
                    if sp.issparse(adata_subset.X):
                        adata_subset.X = adata_subset.X.expm1()  # sparse-safe exp(x)-1
                        # Add small pseudocount to avoid zeros
                        adata_subset.X.data = np.maximum(adata_subset.X.data, 1e-10)
                    else:
                        adata_subset.X = np.expm1(adata_subset.X)
                        # Add small pseudocount to avoid zeros
                        adata_subset.X = np.maximum(adata_subset.X, 1e-10)

                # Normalize to CPM (1e6)
                import scanpy as sc
                sc.pp.normalize_total(adata_subset, target_sum=1e6)
                sc.pp.log1p(adata_subset)

            # Add domain labels for adversarial training
            adata_train_subset.obs['domain'] = 'train'
            adata_test_subset.obs['domain'] = 'test'

            # Combine train and test data for scnym (required for domain adaptation)
            adata_combined = ad.concat([adata_train_subset, adata_test_subset],
                                     join='outer', index_unique='_')

            # Save combined AnnData object
            adata_combined.write_h5ad(combined_h5ad_path)

            # Create Python script for scnym execution
            scnym_script = f'''
import sys
import pickle
import time
import torch
from memory_profiler import memory_usage

try:
    from scnym.api import scnym_api
    import anndata as ad

    # Load combined data (NOT measured)
    adata = ad.read_h5ad('{combined_h5ad_path}')

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        scnym_api(adata=adata, task='train', groupby='Ground_Truth_Celltype',
                  out_path='{temp_folder}', config='no_new_identity')
        scnym_api(adata=adata, task='predict', key_added='scNym',
                  trained_model='{temp_folder}', out_path='{temp_folder}', config='no_new_identity')

        test_mask = adata.obs['domain'] == 'test'
        adata_test_results = adata[test_mask, :].copy()

        if 'scNym' in adata_test_results.obs.columns:
            predictions = list(adata_test_results.obs['scNym'])
        else:
            raise ValueError("scNym prediction column not found in output")
        confidence_scores = list(adata_test_results.obs['scNym_confidence']) if 'scNym_confidence' in adata_test_results.obs.columns else [0.5] * len(predictions)
        true_labels = list(adata_test_results.obs['Ground_Truth_Celltype'])
        cell_ids = list(adata_test_results.obs.index)

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return predictions, confidence_scores, true_labels, cell_ids, _peak_vram_mb

    mem_list, (predictions, confidence_scores, true_labels, cell_ids, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': predictions,
        'confidence_scores': confidence_scores,
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

    print(f"Predictions: {{len(predictions)}} cells, unique: {{set(predictions)}}")

except Exception as e:
    print(f"scnym execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scnym in conda environment
            # Set environment to avoid Jupyter backend conflicts
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend
            env['CUDA_VISIBLE_DEVICES'] = ''  # Force CPU due to RTX 5060 sm_120 incompatibility with PyTorch

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scnym_env",
                    "python",
                    "-c", scnym_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scnym execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scnym results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidence_scores = results['confidence_scores']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Process confidence scores (scnym provides native scores)
            confidence_scores = [
                0.0 if pred == "Unknown" else (max(0.0, min(1.0, float(conf))) if pd.notna(conf) and isinstance(conf, (int, float)) else 0.5)
                for pred, conf in zip(cell_predictions, confidence_scores)
            ]

            print(f"scnym completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")

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
            warnings.warn("scnym execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scnym error: {str(e)}")
            return default_return()


# For backward compatibility
run_scnym = run_scnym_function