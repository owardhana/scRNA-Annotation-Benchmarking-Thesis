# run_scDeepSort.py
#################################################
# scDeepSort Function for Python Benchmarking Framework
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
from typing import Dict, List, Any


def run_scDeepSort_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scDeepSort Cell Type Annotation Function

    Purpose: Run scDeepSort algorithm using pre-trained deep learning models
    Inputs:
      - adata_train: Training AnnData object (used as reference)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scDeepSort's deep learning classifier (https://scdeepsort.readthedocs.io/)
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
    conda_env_path = "/home/oliver/miniconda3/envs/scDeepSort_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scDeepSort_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            train_csv_path = os.path.join(temp_dir, "train_data.csv")
            train_celltypes_path = os.path.join(temp_dir, "train_celltypes.csv")
            test_csv_path = os.path.join(temp_dir, "test_data.csv")
            model_save_path = os.path.join(temp_dir, "model_save")
            results_path = os.path.join(temp_dir, "results")
            results_csv_path = os.path.join(results_path, "predicted_celltype.csv")

            # Create directories
            os.makedirs(model_save_path, exist_ok=True)
            os.makedirs(results_path, exist_ok=True)

            # Get common genes between train and test
            common_genes = sorted(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset to common genes using AnnData native indexing
            adata_train_subset = adata_train[:, common_genes].copy()
            adata_test_subset = adata_test[:, common_genes].copy()

            # Convert to DataFrame (cells x genes) and transpose to (genes x cells) for scDeepSort
            # scDeepSort expects: rows=genes, columns=cells
            train_df = adata_train_subset.to_df().T
            test_df = adata_test_subset.to_df().T

            # Save expression data as CSV
            train_df.to_csv(train_csv_path)
            test_df.to_csv(test_csv_path)

            # Prepare cell types CSV (required format: Cell, Cell_type columns)
            train_celltypes_df = pd.DataFrame({
                'Cell': adata_train_subset.obs.index,
                'Cell_type': adata_train_subset.obs['Ground_Truth_Celltype']
            })
            train_celltypes_df.to_csv(train_celltypes_path, index=True)

            # Create minimal Python script following official API documentation
            scdeepsort_script = f'''
import sys
import os
import pickle
import time
try:
    import torch as _torch
    _has_torch = True
except ImportError:
    _has_torch = False
from memory_profiler import memory_usage

try:
    from deepsort import DeepSortClassifier

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if _has_torch and _torch.cuda.is_available():
            _torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        model = DeepSortClassifier(
            species='human', tissue='Blood', dense_dim=400, hidden_dim=200,
            gpu_id=-1, n_layers=1, random_seed=1, n_epochs=50,
            num_neighbors=100, dropout=0.3, learning_rate=0.0001)

        train_files = [("{train_csv_path}", "{train_celltypes_path}")]
        model.fit(train_files, save_path="{model_save_path}")
        model.predict("{test_csv_path}", save_path="{results_path}", model_path="{model_save_path}")

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = _torch.cuda.max_memory_allocated() if (_has_torch and _torch.cuda.is_available()) else 0
        return _peak_vram_bytes / (1024 * 1024)

    mem_list, peak_vram_mb = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    metadata = {{
        'peak_system_memory_mb': peak_system_memory_mb,
        'peak_vram_mb': peak_vram_mb,
        'method_walltime': method_walltime,
    }}
    with open("{temp_dir}/metadata.pkl", 'wb') as f:
        pickle.dump(metadata, f)

    print(f"Peak system memory: {{peak_system_memory_mb:.2f}} MB, VRAM: {{peak_vram_mb:.2f}} MB, walltime: {{method_walltime:.2f}}s")
    print("scDeepSort execution completed successfully")

except Exception as e:
    print(f"scDeepSort execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scDeepSort in conda environment
            # Set DGL backend environment variable (prevents interactive prompt)
            env = os.environ.copy()
            env['DGLBACKEND'] = 'pytorch'
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend to avoid Jupyter conflicts

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scDeepSort_env",
                    "python",
                    "-c", scdeepsort_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout
            )

            if result.returncode != 0:
                warnings.warn(f"scDeepSort execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # scDeepSort saves results as CSV in the results directory
            # Look for the prediction file
            if not os.path.exists(results_csv_path):
                # Try to find any CSV file in results directory
                result_files = [f for f in os.listdir(results_path) if f.endswith('.csv')]
                if result_files:
                    results_csv_path = os.path.join(results_path, result_files[0])
                else:
                    warnings.warn("Could not find scDeepSort prediction results")
                    return default_return()

            # Read predictions from CSV
            predictions_df = pd.read_csv(results_csv_path, index_col=0)

            # Extract predictions (handle different possible column names)
            if 'predicted_cell_type' in predictions_df.columns:
                cell_predictions = predictions_df['predicted_cell_type'].values
            elif 'prediction' in predictions_df.columns:
                cell_predictions = predictions_df['prediction'].values
            else:
                # Use first column if column name is unclear
                cell_predictions = predictions_df.iloc[:, 0].values

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores (scDeepSort doesn't provide explicit confidence)
            confidence_scores = [0.8 if pred != "Unknown" else 0.0 for pred in cell_predictions]

            # Load timing/memory from metadata file
            metadata_path = os.path.join(temp_dir, "metadata.pkl")
            metadata = {}
            if os.path.exists(metadata_path):
                with open(metadata_path, 'rb') as f:
                    metadata = pickle.load(f)

            print(f"scDeepSort completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
                'confidence_scores': confidence_scores,
                'cell_ids': list(adata_test.obs.index),
                'peak_memory_mb':  metadata.get('peak_system_memory_mb', None),
                'peak_vram_mb':    metadata.get('peak_vram_mb', None),
                'method_walltime': metadata.get('method_walltime', None),
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scDeepSort execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scDeepSort error: {str(e)}")
            return default_return()


# For backward compatibility
run_scDeepSort = run_scDeepSort_function