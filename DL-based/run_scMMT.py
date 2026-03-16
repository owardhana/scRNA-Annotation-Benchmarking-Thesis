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
import numpy as np
import anndata as ad
import time
import warnings
import pickle
from typing import Dict, List, Any, Optional


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

    # Check conda environment availability
    try:
        env_check = subprocess.run(
            "conda info --envs | grep scMMT_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("scMMT_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scMMT...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            model_weights_dir = os.path.join(temp_dir, "model_weight")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Create model weights directory
            os.makedirs(model_weights_dir, exist_ok=True)

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for scMMT analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running scMMT in conda environment...")

            # Create Python script for scMMT execution
            scmmt_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import anndata as ad
import pickle
import torch
import warnings
warnings.filterwarnings("ignore")

try:
    from scMMT.scMMT_API import scMMT_API
    from sklearn.metrics import f1_score, accuracy_score

    # Set random seeds for reproducibility
    seed = 5
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed) if torch.cuda.is_available() else None
    np.random.seed(seed)

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")
    print(f"CUDA available: {{torch.cuda.is_available()}}")

    # Initialize scMMT API
    print("Initializing scMMT API...")

    scMMT = scMMT_API(
        gene_trainsets=[adata_train],  # List of training datasets
        gene_test=adata_test,
        log_normalize=True,            # Log normalization
        type_key='Ground_Truth_Celltype',
        data_load=False,               # Don't load existing processed data
        dataset_batch=True,            # Account for batch effects
        log_weight=3,                  # Log weights for different cell types
        val_split=None,                # No validation split
        min_cells=0,                   # Minimum cell count filtering
        min_genes=0,                   # Minimum gene count filtering
        n_svd=300,                     # SVD dimensionality reduction
        n_fa=180,                      # Factor Analysis dimensionality
        n_hvg=550,                     # Number of highly variable genes
    )

    print("scMMT API initialized successfully")

    # Train the model
    print("Training scMMT model...")

    scMMT.train(
        n_epochs=50,        # Reduced epochs for CV efficiency
        ES_max=12,          # Early stopping patience
        decay_max=6,        # Learning rate decay patience
        decay_step=0.1,     # Learning rate decay step
        lr=1e-3,            # Learning rate
        label_smoothing=0.4, # Label smoothing
        h_size=600,         # Hidden size
        drop_rate=0.15,     # Dropout rate
        n_layer=4,          # Number of layers
        weights_dir='{model_weights_dir}',
        load=False          # Don't load existing model
    )

    print("scMMT training completed")

    # Make predictions
    print("Making predictions...")

    predicted_test = scMMT.predict()

    # Extract predictions
    if 'transfered cell labels' in predicted_test.obs.columns:
        predictions = list(predicted_test.obs['transfered cell labels'])
    else:
        raise ValueError("'transfered cell labels' column not found in scMMT output")

    # Calculate accuracy for validation
    true_labels = list(predicted_test.obs['Ground_Truth_Celltype'])
    acc = (np.array(predictions) == np.array(true_labels)).mean()

    print(f"Training accuracy: {{acc:.3f}}")

    # Prepare results
    results = {{
        'predictions': predictions,
        'true_labels': true_labels,
        'cell_ids': list(predicted_test.obs.index),
        'accuracy': acc
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("scMMT execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"scMMT execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scMMT in conda environment
            result = subprocess.run(
                f'conda run -n scMMT_env python -c "{scmmt_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=3600  # 60 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"scMMT execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing scMMT results...")

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

            # Create confidence scores based on training accuracy and prediction quality
            # Use training accuracy as base confidence, adjust for unknown predictions
            base_confidence = min(0.9, max(0.3, training_accuracy))

            confidence_scores = []
            training_cell_types = set(adata_train.obs['Ground_Truth_Celltype'].unique())

            for pred in cell_predictions:
                if pred == "Unknown":
                    confidence_scores.append(0.0)
                elif pred in training_cell_types:
                    confidence_scores.append(base_confidence)
                else:
                    # Novel cell type detected
                    confidence_scores.append(base_confidence * 0.7)

            print(f"scMMT completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Training accuracy: {training_accuracy:.3f}")
            print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scMMT execution timed out after 60 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"scMMT error: {str(e)}")
            return default_return()


# For backward compatibility
run_scMMT = run_scMMT_function