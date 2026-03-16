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

    # Check conda environment availability
    try:
        env_check = subprocess.run(
            "conda info --envs | grep Cell_BLAST_env",
            shell=True,
            capture_output=True,
            text=True
        )
        if env_check.returncode != 0:
            warnings.warn("Cell_BLAST_env conda environment not found. Please create it first.")
            return default_return()
    except Exception as e:
        warnings.warn(f"Could not check conda environments: {str(e)}")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for Cell_BLAST...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for Cell_BLAST analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            print("Running Cell_BLAST in conda environment...")

            # Create Python script for Cell_BLAST execution
            cell_blast_script = f'''
import sys
import pandas as pd
import numpy as np
import anndata as ad
import pickle

try:
    import Cell_BLAST as cb

    # Load data
    print("Loading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")
    print(f"Training cell types: {{list(adata_train.obs['Ground_Truth_Celltype'].unique())}}")

    # Find variable genes
    print("Finding variable genes...")
    try:
        axes = cb.data.find_variable_genes(adata_train, grouping=None)
        var_genes = adata_train.var.query("variable_genes").index
        print(f"Found {{len(var_genes)}} variable genes")
    except Exception as e:
        print(f"Variable gene detection failed: {{e}}")
        # Fallback: use all genes
        var_genes = adata_train.var.index
        print(f"Using all {{len(var_genes)}} genes")

    # Train DIRECTi model (single model as specified)
    print("Training DIRECTi model...")
    model = cb.directi.fit_DIRECTi(
        adata_train,
        genes=var_genes,
        latent_dim=10,
        cat_dim=20,
        random_seed=0
    )
    print("DIRECTi model trained successfully")

    # Build BLAST database
    print("Building BLAST database...")
    blast = cb.blast.BLAST([model], adata_train)
    print("BLAST database created successfully")

    # Query test cells
    print(f"Querying {{len(adata_test.obs)}} test cells...")
    hits = blast.query(adata_test)
    print(f"Query completed, got {{len(hits)}} hits")

    # Reconcile predictions across models and filter by p-value
    print("Filtering hits by p-value...")
    hits = hits.reconcile_models()
    hits = hits.filter(by="pval", cutoff=0.05)
    print(f"After p-value filtering (cutoff=0.05): {{len(hits)}} hits")

    # Extract predictions and confidence scores
    print("Extracting predictions...")
    predictions = []
    confidences = []
    cell_ids = adata_test.obs.index.tolist()

    for cell_idx, cell_id in enumerate(cell_ids):
        try:
            # Get hits for this cell
            cell_hits = hits[cell_idx]

            if len(cell_hits) == 0:
                # No significant hits
                predictions.append("Unknown")
                confidences.append(0.0)
            else:
                # Get top hit
                top_hit_df = cell_hits.to_data_frames()[cell_id]

                if len(top_hit_df) == 0:
                    predictions.append("Unknown")
                    confidences.append(0.0)
                else:
                    # Extract celltype from top hit
                    top_celltype = top_hit_df.iloc[0]['Ground_Truth_Celltype']
                    top_pval = top_hit_df.iloc[0]['empirical_pval']

                    predictions.append(top_celltype)
                    # Convert p-value to confidence: 1 - pval
                    confidences.append(max(0, 1 - top_pval))
        except Exception as e:
            print(f"WARNING: Failed to extract prediction for cell {{cell_id}}: {{e}}")
            predictions.append("Unknown")
            confidences.append(0.0)

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidences': confidences,
        'true_labels': list(adata_test.obs['Ground_Truth_Celltype']),
        'cell_ids': cell_ids
    }}

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("Cell_BLAST execution completed successfully")
    print(f"Predictions: {{len(predictions)}} cells")
    print(f"Assigned: {{sum([p != 'Unknown' for p in predictions])}} cells")
    print(f"Unknown: {{sum([p == 'Unknown' for p in predictions])}} cells")
    print(f"Unique predictions: {{set(predictions)}}")

except Exception as e:
    print(f"Cell_BLAST execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute Cell_BLAST in conda environment
            result = subprocess.run(
                f'conda run -n Cell_BLAST_env python -c "{cell_blast_script}"',
                shell=True,
                capture_output=True,
                text=True,
                timeout=2400  # 40 minute timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"Cell_BLAST execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            print("Parsing Cell_BLAST results...")

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
                'cell_ids': cell_ids
            }

        except subprocess.TimeoutExpired:
            warnings.warn("Cell_BLAST execution timed out after 40 minutes")
            return default_return()
        except Exception as e:
            warnings.warn(f"Cell_BLAST error: {str(e)}")
            import traceback
            traceback.print_exc()
            return default_return()


# For backward compatibility
run_Cell_BLAST = run_Cell_BLAST_function
