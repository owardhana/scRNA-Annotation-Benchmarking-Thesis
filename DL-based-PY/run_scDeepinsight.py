# run_scDeepinsight.py
#################################################
# scDeepinsight Function for Python Benchmarking Framework
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


def run_scDeepinsight_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scDeepinsight Cell Type Annotation Function

    Purpose: Run scDeepinsight algorithm using image-based classification
    Inputs:
      - adata_train: Training AnnData object (used for reference/format compatibility)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scDeepinsight's image-based deep learning classifier (two-stage process)
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
    conda_env_path = "/home/oliver/miniconda3/envs/scDeepInsight_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scDeepInsight_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scDeepinsight...")

            # Prepare file paths
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            barcode_csv_path = os.path.join(temp_dir, "barcode.csv")
            image_npy_path = os.path.join(temp_dir, "image.npy")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for scDeepinsight analysis")

            # Subset test dataset to common genes and ensure consistent ordering
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save test AnnData object
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running scDeepinsight in conda environment...")

            # Create Python script for scDeepinsight execution
            scdeepinsight_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory (before imports, after data conversion)

import sys
import pickle

try:
    from scdeepinsight import pbmc
    import anndata as ad

    # Load test data
    print("Loading test data...")
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")

    # Stage 1: ImageTransform
    # Convert scRNA-seq data to image representation
    print("Stage 1: Converting scRNA-seq data to image representation...")

    pbmc.ImageTransform(
        query_path='{test_h5ad_path}',
        barcode_path='{barcode_csv_path}',
        image_path='{image_npy_path}'
    )

    print("Image transformation completed")

    # Stage 2: Annotate
    # Use pre-trained model for cell type annotation
    print("Stage 2: Running cell type annotation...")

    pred_labels = pbmc.Annotate(
        barcode_path='{barcode_csv_path}',
        image_path='{image_npy_path}',
        batch_size=128
    )

    print("Annotation completed")

    # Prepare results
    results = {{
        'predictions': pred_labels,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': list(adata_test.obs.index)
    }}

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)  # Convert bytes to MB

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("scDeepinsight execution completed successfully")
    print(f"Predictions: {{len(pred_labels)}} cells")
    print(f"Unique predictions: {{set(pred_labels)}}")

except Exception as e:
    print(f"scDeepinsight execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scDeepinsight in conda environment
            # Set environment to avoid Jupyter backend conflicts
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scDeepInsight_env",
                    "python",
                    "-c", scdeepinsight_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout
            )

            if result.returncode != 0:
                warnings.warn(f"scDeepinsight execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing scDeepinsight results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scDeepinsight results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            peak_memory_mb = results.get('peak_memory_mb', None)  # Extract peak memory

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Create confidence scores (scDeepinsight doesn't provide explicit confidence)
            confidence_scores = [0.8 if pred != "Unknown" else 0.0 for pred in cell_predictions]

            print(f"scDeepinsight completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb  # Add peak memory to return dict
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scDeepinsight execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scDeepinsight error: {str(e)}")
            return default_return()


# For backward compatibility
run_scDeepinsight = run_scDeepinsight_function