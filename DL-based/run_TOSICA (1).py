# run_TOSICA.py
#################################################
# TOSICA Function for Python Benchmarking Framework
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


def run_TOSICA_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    TOSICA Cell Type Annotation Function

    Purpose: Run TOSICA algorithm using pathway-based attention mechanism
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses TOSICA's pathway-based attention with gene set masks
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
    conda_env_path = "/home/oliver/miniconda3/envs/TOSICA_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"TOSICA_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Note: TOSICA will create the project folder itself

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            # Create Python script for TOSICA execution
            tosica_script = f'''
import sys
import os
import pickle
import time
import torch
from memory_profiler import memory_usage

# Change to temp directory before running TOSICA
# TOSICA constructs paths as: os.getcwd() + '/' + project
os.chdir('{temp_dir}')

try:
    import TOSICA
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

        TOSICA.train(
            adata_train,
            gmt_path="human_reactome",
            project='tosica_project',
            label_name="Ground_Truth_Celltype",
            epochs=3
        )

        project_path = 'tosica_project'
        all_files = os.listdir(project_path)
        model_files = sorted([f for f in all_files if f.endswith('.pth')])
        if not model_files:
            raise ValueError(f"No trained model weights (.pth files) found in {{project_path}}. Available files: {{all_files}}")
        model_weight_path = os.path.join(project_path, model_files[-1])

        new_adata = TOSICA.pre(
            adata_test,
            model_weight_path=model_weight_path,
            project='tosica_project',
            laten=True
        )

        if 'Prediction' in new_adata.obs.columns:
            predictions = list(new_adata.obs['Prediction'])
        else:
            raise ValueError("Prediction column not found in TOSICA output")
        probabilities = list(new_adata.obs['Probability']) if 'Probability' in new_adata.obs.columns else [0.5] * len(predictions)
        true_labels = list(new_adata.obs['Ground_Truth_Celltype'])
        cell_ids = list(new_adata.obs.index)

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return predictions, probabilities, true_labels, cell_ids, _peak_vram_mb

    mem_list, (predictions, probabilities, true_labels, cell_ids, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': predictions,
        'probabilities': probabilities,
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
    print(f"TOSICA execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute TOSICA in conda environment
            # Set environment to avoid Jupyter backend conflicts
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "TOSICA_env",
                    "python",
                    "-c", tosica_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"TOSICA execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("TOSICA results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            probabilities = results['probabilities']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Use TOSICA probabilities as confidence scores
            confidence_scores = [
                0.0 if pred == "Unknown" else (max(0.0, min(1.0, float(prob))) if pd.notna(prob) and isinstance(prob, (int, float)) else 0.5)
                for pred, prob in zip(cell_predictions, probabilities)
            ]

            print(f"TOSICA completed successfully:")
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
            warnings.warn("TOSICA execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"TOSICA error: {str(e)}")
            return default_return()


# For backward compatibility
run_TOSICA = run_TOSICA_function