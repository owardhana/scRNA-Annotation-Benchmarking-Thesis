# run_SCANVI.py
#################################################
# SCANVI Function for Python Benchmarking Framework
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


def run_SCANVI_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    SCANVI Cell Type Annotation Function

    Purpose: Run SCANVI algorithm using semi-supervised variational inference
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scvi-tools SCANVI for semi-supervised cell type annotation
               with variational autoencoder and online query mapping
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
    conda_env_path = "/home/oliver/miniconda3/envs/SCANVI_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"SCANVI_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            model_dir = os.path.join(temp_dir, "scanvi_model")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Ensure counts layer exists (SCANVI requires raw counts)
            # Check if counts layer exists, if not use .X
            if 'counts' not in adata_train_subset.layers:
                adata_train_subset.layers['counts'] = adata_train_subset.X.copy()
            if 'counts' not in adata_test_subset.layers:
                adata_test_subset.layers['counts'] = adata_test_subset.X.copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            # Create Python script for SCANVI execution
            scanvi_script = f'''
import sys
import os
import pickle
import warnings
warnings.filterwarnings('ignore')

try:
    import scvi
    import scanpy as sc
    import anndata as ad
    import numpy as np
    import time
    import torch
    from memory_profiler import memory_usage

    # Set scvi settings
    scvi.settings.verbosity = 1
    scvi.settings.progress_bar_style = "tqdm"

    # Load data (NOT measured)
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # --- method execution (measured) ---
    _timing = [0.0]

    def _run_method():
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        _t0 = time.perf_counter()

        # PART 1: Train SCVI reference model
        scvi.model.SCVI.setup_anndata(adata_train, layer="counts")
        scvi_ref = scvi.model.SCVI(adata_train, n_layers=2, n_latent=30, gene_likelihood="nb")
        scvi_ref.train(max_epochs=100)

        # PART 2: Initialize and train SCANVI
        scanvi_ref = scvi.model.SCANVI.from_scvi_model(
            scvi_ref, unlabeled_category="Unknown", labels_key="Ground_Truth_Celltype")
        scanvi_ref.train(max_epochs=20, n_samples_per_label=100)

        # PART 3: Save reference model
        scanvi_ref.save('{model_dir}', overwrite=True)

        # PART 4: Prepare and predict on query data
        scvi.model.SCANVI.prepare_query_anndata(adata_test, '{model_dir}')
        scanvi_query = scvi.model.SCANVI.load_query_data(adata_test, '{model_dir}')
        scanvi_query.train(
            max_epochs=100,
            plan_kwargs={{"weight_decay": 0.0}},
            check_val_every_n_epoch=10,
        )

        # PART 5: Get predictions
        predictions = scanvi_query.predict()
        prediction_probs = scanvi_query.predict(soft=True)
        confidence_scores = prediction_probs.max(axis=1).tolist()
        true_labels = list(adata_test.obs['Ground_Truth_Celltype'])
        cell_ids = list(adata_test.obs.index)

        _timing[0] = time.perf_counter() - _t0
        _peak_vram_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        _peak_vram_mb = _peak_vram_bytes / (1024 * 1024)
        return predictions, confidence_scores, true_labels, cell_ids, _peak_vram_mb

    mem_list, (predictions, confidence_scores, true_labels, cell_ids, peak_vram_mb) = memory_usage(
        (_run_method, [], {{}}), interval=0.1, include_children=True, retval=True)
    peak_system_memory_mb = max(mem_list)
    method_walltime = _timing[0]

    results = {{
        'predictions': list(predictions),
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
    print(f"SCANVI execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute SCANVI in conda environment
            # Set environment - enable GPU
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "SCANVI_env",
                    "python",
                    "-c", scanvi_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"SCANVI execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("SCANVI results file not found")
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

            # Process confidence scores
            confidence_scores = [
                0.0 if pred == "Unknown" else (max(0.0, min(1.0, float(conf))) if pd.notna(conf) and isinstance(conf, (int, float)) else 0.5)
                for pred, conf in zip(cell_predictions, confidence_scores)
            ]

            print(f"SCANVI completed successfully:")
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
            warnings.warn("SCANVI execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"SCANVI error: {str(e)}")
            return default_return()


# For backward compatibility
run_SCANVI = run_SCANVI_function
