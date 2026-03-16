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
import numpy as np
import anndata as ad
import time
import warnings
import pickle
from typing import Dict, List, Any, Optional


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

    # Check conda environment availability
    try:
        env_check = subprocess.run(
            "conda info --envs | grep scnym_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("scnym_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scnym...")

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

            print(f"Using {len(common_genes)} common genes for scnym analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Add domain labels for adversarial training
            adata_train_subset.obs['domain'] = 'train'
            adata_test_subset.obs['domain'] = 'test'

            # Combine train and test data for scnym (required for domain adaptation)
            adata_combined = ad.concat([adata_train_subset, adata_test_subset],
                                     join='outer', index_unique='_')

            # Save combined AnnData object
            adata_combined.write_h5ad(combined_h5ad_path)

            print("Running scnym in conda environment...")

            # Create Python script for scnym execution
            scnym_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    from scnym.api import scnym_api

    # Load combined data
    print("Loading combined data...")
    adata = ad.read_h5ad('{combined_h5ad_path}')

    print(f"Combined data: {{len(adata.obs)}} cells x {{len(adata.var)}} genes")
    print(f"Domains: {{list(adata.obs['domain'].unique())}}")
    print(f"Cell types: {{list(adata.obs['Ground_Truth_Celltype'].unique())}}")

    # Stage 1: Train scnym model
    print("Stage 1: Training scnym model...")

    scnym_api(
        adata=adata,
        task='train',
        groupby='Ground_Truth_Celltype',
        out_path='{temp_folder}',
        config='no_new_identity',
    )

    print("scnym training completed")

    # Stage 2: Predict with scnym
    print("Stage 2: Making predictions with scnym...")

    scnym_api(
        adata=adata,
        task='predict',
        key_added='scNym',
        trained_model='{temp_folder}',
        out_path='{temp_folder}',
        config='no_new_identity',
    )

    print("scnym prediction completed")

    # Extract test data predictions
    test_mask = adata.obs['domain'] == 'test'
    adata_test_results = adata[test_mask, :].copy()

    # Extract predictions and confidence scores
    if 'scNym' in adata_test_results.obs.columns:
        predictions = list(adata_test_results.obs['scNym'])
    else:
        raise ValueError("scNym prediction column not found in output")

    if 'scNym_confidence' in adata_test_results.obs.columns:
        confidence_scores = list(adata_test_results.obs['scNym_confidence'])
    else:
        # Use default confidence scores if not available
        confidence_scores = [0.5] * len(predictions)

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidence_scores': confidence_scores,
        'true_labels': list(adata_test_results.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test_results.obs.index)
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("scnym execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"scnym execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scnym in conda environment
            result = subprocess.run(
                f'conda run -n scnym_env python -c "{scnym_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=2400  # 40 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scnym execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing scnym results...")

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

            # Process confidence scores
            processed_confidence_scores = []
            for conf in confidence_scores:
                if pd.notna(conf) and isinstance(conf, (int, float)):
                    # Ensure confidence is between 0 and 1
                    conf_val = max(0.0, min(1.0, float(conf)))
                    processed_confidence_scores.append(conf_val)
                else:
                    processed_confidence_scores.append(0.5)  # Default confidence

            # Adjust confidence for Unknown predictions
            for i, pred in enumerate(cell_predictions):
                if pred == "Unknown":
                    processed_confidence_scores[i] = 0.0

            print(f"scnym completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Confidence range: {min(processed_confidence_scores):.3f} - {max(processed_confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': processed_confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scnym execution timed out after 40 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"scnym error: {str(e)}")
            return default_return()


# For backward compatibility
run_scnym = run_scnym_function