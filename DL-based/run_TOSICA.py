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
import numpy as np
import anndata as ad
import time
import warnings
import pickle
from typing import Dict, List, Any, Optional


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

    # Check conda environment availability
    try:
        env_check = subprocess.run(
            "conda info --envs | grep TOSICA_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("TOSICA_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for TOSICA...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            project_folder = os.path.join(temp_dir, "tosica_project")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for TOSICA analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running TOSICA in conda environment...")

            # Create Python script for TOSICA execution
            tosica_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import TOSICA

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Stage 1: Train TOSICA model
    print("Stage 1: Training TOSICA model...")

    TOSICA.train(
        adata_train,
        gmt_path="human_reactome",
        project='{project_folder}',
        label_name="Ground_Truth_Celltype"
    )

    print("TOSICA training completed")

    # Find the trained model weight file
    model_files = []
    if os.path.exists('{project_folder}'):
        for file in os.listdir('{project_folder}'):
            if file.startswith('model-') and file.endswith('.pth'):
                model_files.append(file)

    if not model_files:
        raise ValueError("No trained model weights found")

    # Use the last model file (highest epoch number)
    model_files.sort()
    model_weight_path = os.path.join('{project_folder}', model_files[-1])
    print(f"Using model weights: {{model_weight_path}}")

    # Stage 2: Predict with TOSICA
    print("Stage 2: Making predictions with TOSICA...")

    new_adata = TOSICA.pre(
        adata_test,
        model_weight_path=model_weight_path,
        project='{project_folder}'
    )

    print("TOSICA prediction completed")

    # Extract predictions and probabilities
    if 'Prediction' in new_adata.obs.columns:
        predictions = list(new_adata.obs['Prediction'])
    else:
        raise ValueError("Prediction column not found in TOSICA output")

    if 'Probability' in new_adata.obs.columns:
        probabilities = list(new_adata.obs['Probability'])
    else:
        # Use default probabilities if not available
        probabilities = [0.5] * len(predictions)

    # Prepare results
    results = {{
        'predictions': predictions,
        'probabilities': probabilities,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index)
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("TOSICA execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"TOSICA execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute TOSICA in conda environment
            result = subprocess.run(
                f'conda run -n TOSICA_env python -c "{tosica_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=2400  # 40 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"TOSICA execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing TOSICA results...")

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

            # Use TOSICA probabilities as confidence scores, normalize if needed
            confidence_scores = []
            for prob in probabilities:
                if pd.notna(prob) and isinstance(prob, (int, float)):
                    # Ensure probability is between 0 and 1
                    conf = max(0.0, min(1.0, float(prob)))
                    confidence_scores.append(conf)
                else:
                    confidence_scores.append(0.5)  # Default confidence

            # Handle predictions that are "Unknown"
            for i, pred in enumerate(cell_predictions):
                if pred == "Unknown":
                    confidence_scores[i] = 0.0

            print(f"TOSICA completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("TOSICA execution timed out after 40 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"TOSICA error: {str(e)}")
            return default_return()


# For backward compatibility
run_TOSICA = run_TOSICA_function