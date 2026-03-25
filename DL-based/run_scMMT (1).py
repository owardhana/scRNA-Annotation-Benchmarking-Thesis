# run_scMMT.py
#################################################
# scMMT Function for Python Benchmarking Framework
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


def run_scMMT_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scMMT Cell Type Annotation Function

    Purpose: Run scMMT algorithm using multi-modal transformer for cross-dataset transfer
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scMMT's multi-modal transformer with cross-dataset adaptation
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
    conda_env_path = "/home/oliver/miniconda3/envs/scMMT_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scMMT_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            model_weights_dir = os.path.join(temp_dir, "model_weight")
            data_dir_path = os.path.join(temp_dir, "scmmt_preprocessed.pkl")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # scMMT works best with raw counts + its own normalization pipeline.
            # Try to recover raw counts so scMMT can normalize properly (normalize_total -> log1p -> scale).
            import scipy.sparse as sp
            use_log_normalize = False
            if 'counts' in adata_train_subset.layers:
                adata_train_subset.X = adata_train_subset.layers['counts'].copy()
                adata_test_subset.X = adata_test_subset.layers['counts'].copy()
                use_log_normalize = True
            elif hasattr(adata_train_subset, 'raw') and adata_train_subset.raw is not None:
                try:
                    raw_train = adata_train_subset.raw.to_adata()
                    raw_test = adata_test_subset.raw.to_adata()
                    raw_common = sorted(set(raw_train.var.index) & set(raw_test.var.index) & set(sorted(common_genes)))
                    if len(raw_common) >= 100:
                        adata_train_subset = raw_train[:, raw_common].copy()
                        adata_test_subset = raw_test[:, raw_common].copy()
                        adata_train_subset.obs['Ground_Truth_Celltype'] = adata_train.obs['Ground_Truth_Celltype'].values
                        adata_test_subset.obs['Ground_Truth_Celltype'] = adata_test.obs['Ground_Truth_Celltype'].values
                        use_log_normalize = True
                except Exception:
                    pass  # Fall back to pre-normalized data
            else:
                # Reverse log1p to approximate raw counts
                if sp.issparse(adata_train_subset.X):
                    adata_train_subset.X = sp.csr_matrix(np.expm1(adata_train_subset.X.toarray()))
                    adata_test_subset.X = sp.csr_matrix(np.expm1(adata_test_subset.X.toarray()))
                else:
                    adata_train_subset.X = np.expm1(np.array(adata_train_subset.X))
                    adata_test_subset.X = np.expm1(np.array(adata_test_subset.X))
                use_log_normalize = True

            # Save AnnData objects for inter-process communication
            adata_train_subset.layers.clear()
            adata_test_subset.layers.clear()
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)
            # Remove /layers group from h5ad — old anndata can't parse encoding-type 'dict'
            for h5_path in [train_h5ad_path, test_h5ad_path]:
                with h5py.File(h5_path, 'a') as f:
                    if 'layers' in f:
                        del f['layers']

            # Create Python script for scMMT execution
            scmmt_script = f'''
import sys
import pickle
import warnings
warnings.filterwarnings("ignore")
import time
import torch
from memory_profiler import memory_usage

try:
    from scMMT.scMMT_API import scMMT_API
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    # Set random seeds for reproducibility
    seed = 5
    torch.manual_seed(seed); torch.cuda.manual_seed_all(seed) if torch.cuda.is_available() else None; np.random.seed(seed)

    # Load data (NOT measured)
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # Data cleaning (NOT measured)
    def clean_data(X):
        if sp.issparse(X):
            X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
            X.data = np.clip(X.data, -1e10, 1e10)
        else:
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            X = np.clip(X, -1e10, 1e10)
        return X

    adata_train.X = clean_data(adata_train.X)
    adata_test.X = clean_data(adata_test.X)

    def has_bad_values(X):
        data = X.data if sp.issparse(X) else X.ravel()
        return np.any(np.isnan(data)) or np.any(np.isinf(data))

    if has_bad_values(adata_train.X) or has_bad_values(adata_test.X):
        raise ValueError("Failed to clean data - still contains NaN or Inf values")

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        scMMT = scMMT_API(
            gene_trainsets=[adata_train],
            gene_test=adata_test,
            log_normalize={use_log_normalize},
            select_hvg=False,
            gene_normalize=True,
            cell_normalize=False,
            type_key='Ground_Truth_Celltype',
            data_load=False,
            data_dir='{data_dir_path}',
            dataset_batch=True,
            log_weight=None,
            val_split=None,
            min_cells=0,
            min_genes=0,
            n_svd=300,
            n_fa=180,
            n_hvg=550,
        )
        scMMT.train(
            n_epochs=200, ES_max=12, decay_max=6, decay_step=0.1,
            lr=1e-3, label_smoothing=0.1, h_size=600, drop_rate=0.15,
            n_layer=4, weights_dir='{model_weights_dir}', load=False)
        predicted_test = scMMT.predict()

        if 'transfered cell labels' in predicted_test.obs.columns:
            predictions = list(predicted_test.obs['transfered cell labels'])
        else:
            raise ValueError("'transfered cell labels' column not found in scMMT output")
        true_labels = list(predicted_test.obs['Ground_Truth_Celltype'])
        cell_ids = list(predicted_test.obs.index)
        acc = (np.array(predictions) == np.array(true_labels)).mean()
        print(f"Training accuracy: {{acc:.3f}}")

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return predictions, true_labels, cell_ids, acc, _peak_vram_mb

    mem_list, (predictions, true_labels, cell_ids, acc, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': predictions,
        'true_labels': true_labels,
        'cell_ids': cell_ids,
        'accuracy': acc,
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
    print(f"scMMT execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scMMT in conda environment
            # Set environment to use non-interactive matplotlib backend
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scMMT_env",
                    "python",
                    "-c", scmmt_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scMMT execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scMMT results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            training_accuracy = results.get('accuracy', 0.0)

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores based on training accuracy
            base_confidence = min(0.9, max(0.3, training_accuracy))
            confidence_scores = [base_confidence if pred != "Unknown" else 0.0 for pred in cell_predictions]

            print(f"scMMT completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")
            print(f"  - Training accuracy: {training_accuracy:.3f}")

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
            warnings.warn("scMMT execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scMMT error: {str(e)}")
            return default_return()


# For backward compatibility
run_scMMT = run_scMMT_function