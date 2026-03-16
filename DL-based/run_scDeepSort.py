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
import time
import warnings
from typing import Dict, List, Any, Optional


def run_scDeepSort_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scDeepSort Cell Type Annotation Function

    Purpose: Run scDeepSort algorithm using pre-trained deep learning models
    Inputs:
      - adata_train: Training AnnData object (used as reference for format compatibility)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scDeepSort's deep learning classifier with conda environment activation
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
            "conda info --envs | grep scDeepSort_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("scDeepSort_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scDeepSort...")

            # Prepare file paths
            train_csv_path = os.path.join(temp_dir, "adata_train.csv")
            train_celltypes_path = os.path.join(temp_dir, "adata_train_celltypes.csv")
            test_csv_path = os.path.join(temp_dir, "adata_test.csv")
            model_save_path = os.path.join(temp_dir, "model_save")
            results_path = os.path.join(temp_dir, "results")

            # Create model and results directories
            os.makedirs(model_save_path, exist_ok=True)
            os.makedirs(results_path, exist_ok=True)

            # Get expression matrices (use raw counts if available, otherwise use X)
            if adata_train.raw is not None:
                train_expr = adata_train.raw.X
                train_genes = adata_train.raw.var.index
            else:
                train_expr = adata_train.X
                train_genes = adata_train.var.index

            if adata_test.raw is not None:
                test_expr = adata_test.raw.X
                test_genes = adata_test.raw.var.index
            else:
                test_expr = adata_test.X
                test_genes = adata_test.var.index

            # Convert sparse matrices to dense if needed
            if hasattr(train_expr, 'todense'):
                train_expr = train_expr.todense()
            if hasattr(test_expr, 'todense'):
                test_expr = test_expr.todense()

            # Convert to numpy arrays
            train_expr = np.array(train_expr)
            test_expr = np.array(test_expr)

            # Get common genes between train and test
            common_genes = list(set(train_genes) & set(test_genes))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            # Subset to common genes
            train_gene_idx = [i for i, gene in enumerate(train_genes) if gene in common_genes]
            test_gene_idx = [i for i, gene in enumerate(test_genes) if gene in common_genes]

            train_expr_subset = train_expr[train_gene_idx, :]
            test_expr_subset = test_expr[test_gene_idx, :]

            # Create gene order list for consistent ordering
            common_genes_ordered = [train_genes[i] for i in train_gene_idx]

            print(f"Using {len(common_genes_ordered)} common genes for scDeepSort analysis")

            # Prepare training data CSV (genes x cells format as required by scDeepSort)
            train_df = pd.DataFrame(
                train_expr_subset,
                index=common_genes_ordered,
                columns=adata_train.obs.index
            )
            train_df.to_csv(train_csv_path)

            # Prepare training cell types CSV
            train_celltypes_df = pd.DataFrame({
                '': range(1, len(adata_train.obs) + 1),
                'Cell': adata_train.obs.index,
                'Cell_type': adata_train.obs['Ground_Truth_Celltype']
            })
            train_celltypes_df.to_csv(train_celltypes_path, index=False)

            # Prepare test data CSV (genes x cells format)
            test_df = pd.DataFrame(
                test_expr_subset,
                index=common_genes_ordered,
                columns=adata_test.obs.index
            )
            test_df.to_csv(test_csv_path)

            print("Running scDeepSort in conda environment...")

            # Create Python script for scDeepSort execution
            scdeepsort_script = f'''
import sys
import os
sys.path.append(os.getcwd())

try:
    from deepsort import DeepSortClassifier

    # Define the model
    model = DeepSortClassifier(
        species='human',
        tissue='Blood',
        dense_dim=50,
        hidden_dim=20,
        gpu_id=0,
        n_layers=2,
        random_seed=1,
        n_epochs=20
    )

    # Prepare file lists
    train_files = [('{train_csv_path}', '{train_celltypes_path}')]
    test_files = ['{test_csv_path}']

    # Fit the model
    print("Training scDeepSort model...")
    model.fit(train_files, save_path='{model_save_path}')

    # Use the saved model to predict
    print("Making predictions...")
    model.predict(test_files, save_path='{results_path}', model_path='{model_save_path}')

    print("scDeepSort execution completed successfully")

except Exception as e:
    print(f"scDeepSort execution failed: {{str(e)}}")
    sys.exit(1)
'''

            # Execute scDeepSort in conda environment
            result = subprocess.run(
                f'conda run -n scDeepSort_env python -c "{scdeepsort_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minute timeout
            )

            if result.returncode != 0:
                warnings.warn(f"scDeepSort execution failed: {result.stderr}")
                return default_return()

            print("Parsing scDeepSort results...")

            # Parse results - scDeepSort typically saves results as CSV files
            # Look for prediction files in results directory
            result_files = os.listdir(results_path)
            prediction_file = None

            for file in result_files:
                if 'predict' in file.lower() or 'result' in file.lower():
                    prediction_file = os.path.join(results_path, file)
                    break

            if prediction_file is None or not os.path.exists(prediction_file):
                # Try to find any CSV file in results
                csv_files = [f for f in result_files if f.endswith('.csv')]
                if csv_files:
                    prediction_file = os.path.join(results_path, csv_files[0])
                else:
                    warnings.warn("Could not find scDeepSort prediction results")
                    return default_return()

            # Read predictions
            try:
                predictions_df = pd.read_csv(prediction_file, index_col=0)

                # Extract predictions - format may vary
                if 'predicted_cell_type' in predictions_df.columns:
                    cell_predictions = predictions_df['predicted_cell_type'].values
                elif 'prediction' in predictions_df.columns:
                    cell_predictions = predictions_df['prediction'].values
                elif len(predictions_df.columns) == 1:
                    cell_predictions = predictions_df.iloc[:, 0].values
                else:
                    # Assume first column contains predictions
                    cell_predictions = predictions_df.iloc[:, 0].values

            except Exception as e:
                warnings.warn(f"Could not parse scDeepSort results: {str(e)}")
                return default_return()

            # Get true labels for test data
            true_labels = list(adata_test.obs['Ground_Truth_Celltype'])

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores (scDeepSort doesn't provide explicit confidence)
            # Use heuristic: high confidence for known predictions, low for Unknown
            confidence_scores = [0.8 if pred != "Unknown" else 0.0 for pred in cell_predictions]

            print(f"scDeepSort completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': list(adata_test.obs.index)
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scDeepSort execution timed out after 30 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"scDeepSort error: {str(e)}")
            return default_return()


# For backward compatibility
run_scDeepSort = run_scDeepSort_function