# run_CIForm.py
#################################################
# CIForm Function for Python Benchmarking Framework
# Input: Train/test AnnData objects and markers from rank_genes_groups
# Output: Standardized results format for CV framework
#################################################

import subprocess
import tempfile
import os
import shutil
import pandas as pd
import anndata as ad
import warnings
import pickle
from typing import Dict, List, Any


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

    # Check conda environment availability (direct directory check)
    conda_env_path = "/home/oliver/miniconda3/envs/CIForm_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"CIForm_env conda environment not found at {conda_env_path}. Please create it first.")
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

            # Validate no cells were lost during subsetting
            if len(adata_train_subset.obs) != len(adata_train.obs):
                warnings.warn(f"Train cells lost: {len(adata_train.obs)} → {len(adata_train_subset.obs)}")
            if len(adata_test_subset.obs) != len(adata_test.obs):
                warnings.warn(f"Test cells lost: {len(adata_test.obs)} → {len(adata_test_subset.obs)}")

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            # Copy helper file to temp directory for subprocess access
            helper_src = os.path.join(os.path.dirname(__file__), "run_CIForm_helper.py")
            helper_dst = os.path.join(temp_dir, "run_CIForm_helper.py")
            shutil.copy(helper_src, helper_dst)

            # Create Python script for CIForm execution using the helper
            ciform_script = f'''
import sys
import pickle
import time
import torch
from memory_profiler import memory_usage

try:
    from run_CIForm_helper import ciForm
    import anndata as ad

    # Load data (NOT measured)
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        predictions = ciForm(
            s=1024,
            Train_adata=adata_train,
            train_labels=adata_train.obs['Ground_Truth_Celltype'].tolist(),
            Test_adata=adata_test,
            n_epochs=50
        )

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return list(predictions), list(adata_test.obs['Ground_Truth_Celltype']), list(adata_test.obs.index), _peak_vram_mb

    mem_list, (predictions, true_labels, cell_ids, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': predictions,
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
    print(f"CIForm execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute CIForm in conda environment with PYTHONPATH set
            env = os.environ.copy()
            env['PYTHONPATH'] = temp_dir
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend to avoid Jupyter conflicts

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "CIForm_env",
                    "python",
                    "-c", ciform_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout
            )

            if result.returncode != 0:
                warnings.warn(f"CIForm execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

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

            # Convert both predictions and true labels to lowercase for case-insensitive comparison
            # CIForm returns lowercase predictions, so we need to normalize both for accurate metrics
            cell_predictions = [pred.lower() for pred in cell_predictions]
            true_labels = [str(label).lower() for label in true_labels]

            # Create confidence scores (CIForm doesn't provide explicit confidence)
            confidence_scores = [0.7 if pred != "unknown" else 0.0 for pred in cell_predictions]

            print(f"CIForm completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'unknown')}")
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
            warnings.warn("CIForm execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"CIForm error: {str(e)}")
            return default_return()


# For backward compatibility
run_CIForm = run_CIForm_function