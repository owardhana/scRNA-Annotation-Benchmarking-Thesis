# run_Cell_BLAST.py
#################################################
# Cell_BLAST Function for Python Benchmarking Framework
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


def run_Cell_BLAST_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    Cell_BLAST Cell Type Annotation Function

    Purpose: Run Cell_BLAST algorithm using DIRECTi model training and BLAST querying
    Inputs:
      - adata_train: Training AnnData object (used for DIRECTi model training)
      - adata_test: Test AnnData object to predict via BLAST query
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Trains DIRECTi model, creates BLAST database, queries test cells
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
    conda_env_path = "/home/oliver/miniconda3/envs/Cell_BLAST_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"Cell_BLAST_env conda environment not found at {conda_env_path}. Please create it first.")
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

            # Cell_BLAST requires raw count data (non-negative values)
            # Check if raw counts are available, otherwise ensure non-negativity

            # Try to get raw counts
            if hasattr(adata_train_subset, 'raw') and adata_train_subset.raw is not None:
                # Extract raw data for common genes
                train_raw = adata_train_subset.raw.to_adata()
                test_raw = adata_test_subset.raw.to_adata()
                # Subset to common genes
                adata_train_subset.X = train_raw[:, sorted(common_genes)].X
                adata_test_subset.X = test_raw[:, sorted(common_genes)].X
            elif 'counts' in adata_train_subset.layers:
                adata_train_subset.X = adata_train_subset.layers['counts']
                adata_test_subset.X = adata_test_subset.layers['counts']
            else:
                # Check if data has negative values (likely normalized)
                if hasattr(adata_train_subset.X, 'min'):
                    min_val = adata_train_subset.X.min()
                else:
                    min_val = adata_train_subset.X.toarray().min() if hasattr(adata_train_subset.X, 'toarray') else np.min(adata_train_subset.X)

                if min_val < 0:
                    warnings.warn("Data contains negative values. Cell_BLAST works best with raw counts. "
                                "Applying expm1 to reverse log1p normalization.")
                    # Try to reverse log1p transformation
                    import scipy.sparse as sp
                    if hasattr(adata_train_subset.X, 'toarray'):
                        adata_train_subset.X = sp.csr_matrix(np.expm1(adata_train_subset.X.toarray()))
                        adata_test_subset.X = sp.csr_matrix(np.expm1(adata_test_subset.X.toarray()))
                    else:
                        adata_train_subset.X = np.expm1(adata_train_subset.X)
                        adata_test_subset.X = np.expm1(adata_test_subset.X)

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running Cell_BLAST in conda environment...")

            # Create Python script for Cell_BLAST execution
            cell_blast_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory (before imports, after data conversion)

import sys
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import Cell_BLAST as cb

    # Load data
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    # Find variable genes
    try:
        axes = cb.data.find_variable_genes(adata_train, grouping=None)
        var_genes = adata_train.var.query("variable_genes").index
    except Exception as e:
        # Fallback: use all genes
        var_genes = adata_train.var.index

    # Train DIRECTi model (single model as specified)
    model = cb.directi.fit_DIRECTi(
        adata_train,
        genes=var_genes,
        latent_dim=10,
        cat_dim=20,
        random_seed=0
    )
    # Build BLAST database and query
    blast = cb.blast.BLAST([model], adata_train)
    hits = blast.query(adata_test)

    # Annotate using majority voting on hits

    # Annotate using the Ground_Truth_Celltype field
    # min_hits=1 means we accept cells with at least 1 hit
    # majority_threshold=0.5 means we need >50% agreement for confident annotation
    annotation_df = hits.annotate(
        field='Ground_Truth_Celltype',
        min_hits=1,
        majority_threshold=0.5,
        return_evidence=True
    )

    # Extract predictions from annotation DataFrame
    # The annotation_df should have a 'Ground_Truth_Celltype' column with predictions
    cell_ids = adata_test.obs.index.tolist()
    predictions = []
    confidences = []

    for cell_id in cell_ids:
        try:
            if cell_id in annotation_df.index:
                pred = annotation_df.loc[cell_id, 'Ground_Truth_Celltype']

                # Handle rejected/ambiguous annotations
                if pred in ['rejected', 'ambiguous']:
                    predictions.append("Unknown")
                    confidences.append(0.0)
                else:
                    predictions.append(str(pred))

                    # Use majority_frac as confidence if available
                    if 'majority_frac' in annotation_df.columns:
                        maj_frac = annotation_df.loc[cell_id, 'majority_frac']
                        confidences.append(float(maj_frac) if not pd.isna(maj_frac) else 0.5)
                    else:
                        # Default confidence
                        confidences.append(0.8)
            else:
                predictions.append("Unknown")
                confidences.append(0.0)
        except Exception as e:
            print(f"WARNING: Failed to extract prediction for cell {{cell_id}}: {{e}}")
            predictions.append("Unknown")
            confidences.append(0.0)

    assigned = sum(1 for p in predictions if p != "Unknown")
    print(f"Predictions: {{assigned}}/{{len(predictions)}} cells assigned ({{len(predictions) - assigned}} Unknown)")

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidences': confidences,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': cell_ids
    }}

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)  # Convert bytes to MB

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print(f"Predictions: {{len(predictions)}} cells, assigned: {{sum([p != 'Unknown' for p in predictions])}}, unique: {{set(predictions)}}")

except Exception as e:
    print(f"Cell_BLAST execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute Cell_BLAST in conda environment
            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "Cell_BLAST_env",
                    "python",
                    "-c", cell_blast_script
                ],
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"Cell_BLAST execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("Cell_BLAST results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidences = results['confidences']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            peak_memory_mb = results.get('peak_memory_mb', None)  # Extract peak memory

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Ensure confidence scores are valid floats
            confidence_scores = []
            for conf in confidences:
                if pd.notna(conf) and isinstance(conf, (int, float)):
                    # Ensure confidence is between 0 and 1
                    confidence_scores.append(max(0.0, min(1.0, float(conf))))
                else:
                    confidence_scores.append(0.0)

            # Handle predictions that are "Unknown"
            for i, pred in enumerate(cell_predictions):
                if pred == "Unknown":
                    confidence_scores[i] = 0.0

            print(f"Cell_BLAST completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Assigned predictions: {sum(1 for p in cell_predictions if p != 'Unknown')}")
            print(f"  - Unique predictions: {set(cell_predictions)}")
            print(f"  - Confidence range: {min(confidence_scores):.3f} - {max(confidence_scores):.3f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb  # Add peak memory to return dict
            }

        except subprocess.TimeoutExpired:
            warnings.warn("Cell_BLAST execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"Cell_BLAST error: {str(e)}")
            import traceback
            traceback.print_exc()
            return default_return()


# For backward compatibility
run_Cell_BLAST = run_Cell_BLAST_function
