#!/usr/bin/env python3
# run_scMatch.py
#################################################
# Benchmarking wrapper for scMatch
# Input: Liu Dataset (70% train, 30% test split)
# Output: Standardized results format with comprehensive metrics
#################################################

import pandas as pd
import numpy as np
import time
import os
import subprocess
import sys
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Activate pandas conversion
pandas2ri.activate()

# Import R packages for data loading
base = importr('base')
seurat = importr('Seurat')

def calculate_metrics(predicted, true_labels):
    """Calculate comprehensive performance metrics"""
    # Overall accuracy
    overall_accuracy = accuracy_score(true_labels, predicted)
    
    # Per-class metrics
    unique_labels = np.unique(np.concatenate([predicted, true_labels]))
    unique_labels = unique_labels[(unique_labels != 'Unknown') & (unique_labels != 'unassigned')]
    
    # Calculate precision, recall, F1 for each class
    precision, recall, f1, support = precision_recall_fscore_support(
        true_labels, predicted, labels=unique_labels, average=None, zero_division=0
    )
    
    # Macro F1 (unweighted average)
    macro_f1 = np.mean(f1)
    
    per_class_metrics = pd.DataFrame({
        'cell_type': unique_labels,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'support': support
    })
    
    return {
        'overall_accuracy': overall_accuracy,
        'macro_f1': macro_f1,
        'per_class_metrics': per_class_metrics
    }

def load_seurat_object(seurat_path):
    """Load Seurat object from R and convert to pandas DataFrame"""
    print("Loading Seurat object...")
    
    # Load Seurat object in R
    ro.r(f'seurat_obj <- readRDS("{seurat_path}")')
    
    # Extract normalized expression data
    ro.r('expr_data <- GetAssayData(seurat_obj, layer = "data")')
    ro.r('expr_data <- as.matrix(expr_data)')
    
    # Get expression matrix as pandas DataFrame  
    expr_matrix = ro.r('expr_data')
    expr_df = pd.DataFrame(np.array(expr_matrix), 
                          index=ro.r('rownames(expr_data)'), 
                          columns=ro.r('colnames(expr_data)'))
    
    # Extract metadata
    ro.r('metadata <- seurat_obj@meta.data')
    metadata = ro.r('metadata')
    metadata_df = pandas2ri.rpy2py(metadata)
    
    # Check for Ground_Truth_Celltype
    if 'Ground_Truth_Celltype' not in metadata_df.columns:
        raise ValueError("Ground_Truth_Celltype not found in Seurat object metadata!")
    
    return expr_df, metadata_df

def stratified_split(expr_df, metadata_df, test_size=0.3, random_state=42):
    """Perform stratified 70/30 split"""
    print("Performing 70/30 stratified split...")
    
    # Get cell type labels
    cell_types = metadata_df['Ground_Truth_Celltype']
    
    # Perform stratified split
    train_indices, test_indices = train_test_split(
        range(len(cell_types)), 
        test_size=test_size, 
        stratify=cell_types, 
        random_state=random_state
    )
    
    # Split expression data
    train_expr = expr_df.iloc[:, train_indices]
    test_expr = expr_df.iloc[:, test_indices]
    
    # Split metadata
    train_meta = metadata_df.iloc[train_indices]
    test_meta = metadata_df.iloc[test_indices]
    
    print(f"Training set: {train_expr.shape[1]} cells")
    print(f"Test set: {test_expr.shape[1]} cells")
    print("Training set cell type distribution:")
    print(train_meta['Ground_Truth_Celltype'].value_counts())
    print("Test set cell type distribution:")
    print(test_meta['Ground_Truth_Celltype'].value_counts())
    
    return train_expr, test_expr, train_meta, test_meta

def prepare_scmatch_data(expr_df, metadata_df, output_dir):
    """Prepare data in scMatch format"""
    print("Preparing data for scMatch...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # scMatch expects genes as rows, cells as columns (transpose not needed)
    # Add cell type information as the last column for reference data
    reference_data = expr_df.T.copy()  # Transpose: cells as rows, genes as columns
    reference_data['cell_type'] = metadata_df['Ground_Truth_Celltype'].values
    
    # Save reference data
    ref_file = os.path.join(output_dir, "reference_data.csv")
    reference_data.to_csv(ref_file, index=True)
    
    return ref_file

def create_custom_reference(train_expr, train_meta, output_dir):
    """Create custom reference database from training data"""
    print("Creating custom reference database...")
    
    # Create reference directory structure
    ref_dir = os.path.join(output_dir, "custom_reference")
    os.makedirs(ref_dir, exist_ok=True)
    
    # Prepare reference data: genes as rows, cells as columns
    # This is the format scMatch expects
    ref_data = train_expr.copy()  # Already in correct format
    
    # Save expression data
    expr_file = os.path.join(ref_dir, "expression.csv")
    ref_data.to_csv(expr_file, index=True)
    
    # Create cell type annotation file
    cell_types = train_meta['Ground_Truth_Celltype']
    annotation_file = os.path.join(ref_dir, "annotation.csv")
    
    # scMatch format: cell_id, cell_type
    annotation_df = pd.DataFrame({
        'cell_id': cell_types.index,
        'cell_type': cell_types.values
    })
    annotation_df.to_csv(annotation_file, index=False)
    
    return ref_dir

def run_scmatch(query_data, reference_dir, output_dir):
    """Run scMatch using subprocess"""
    print("Running scMatch...")
    
    # Save query data
    query_file = os.path.join(output_dir, "query_data.csv")
    query_data.to_csv(query_file, index=True)
    
    try:
        # Basic scMatch command - using custom reference
        # Note: This assumes scMatch.py is available in PATH or current directory
        # You may need to modify the path to scMatch.py
        cmd = [
            "python", "scMatch.py",
            "--refDS", reference_dir,
            "--dFormat", "csv",
            "--testDS", query_file,
            "--coreNum", "4"
        ]
        
        print(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            cwd=output_dir,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode != 0:
            print(f"scMatch failed with error: {result.stderr}")
            return None
        
        # Look for output file - scMatch typically creates Results_*.xlsx
        output_files = [f for f in os.listdir(output_dir) if f.startswith("Results_") and f.endswith(".xlsx")]
        
        if output_files:
            return os.path.join(output_dir, output_files[0])
        else:
            print("No scMatch output file found")
            return None
            
    except subprocess.TimeoutExpired:
        print("scMatch timed out")
        return None
    except FileNotFoundError:
        print("scMatch.py not found. Please ensure scMatch is installed and accessible.")
        return None
    except Exception as e:
        print(f"Error running scMatch: {str(e)}")
        return None

def parse_scmatch_results(output_file, test_meta):
    """Parse scMatch Excel output"""
    print("Parsing scMatch results...")
    
    try:
        # Read Excel output
        results_df = pd.read_excel(output_file)
        
        # scMatch output format may vary, adapt as needed
        # Typically contains cell IDs and predicted cell types
        if 'Predicted_Cell_Type' in results_df.columns:
            predicted_col = 'Predicted_Cell_Type'
        elif 'predicted_cell_type' in results_df.columns:
            predicted_col = 'predicted_cell_type'
        else:
            # Try to find any column that might contain predictions
            potential_cols = [col for col in results_df.columns if 'predict' in col.lower() or 'type' in col.lower()]
            if potential_cols:
                predicted_col = potential_cols[0]
            else:
                print("Could not find prediction column in scMatch output")
                return None
        
        # Create standardized output
        predictions = results_df[predicted_col].fillna('Unknown')
        
        # Map to test cells (may need adjustment based on scMatch output format)
        if 'Cell_ID' in results_df.columns:
            cell_ids = results_df['Cell_ID']
        else:
            cell_ids = results_df.index
            
        return predictions, cell_ids
        
    except Exception as e:
        print(f"Error parsing scMatch results: {str(e)}")
        return None, None

def create_dummy_results(test_meta):
    """Create dummy results if scMatch fails"""
    print("Creating dummy results (scMatch failed)...")
    
    n_cells = len(test_meta)
    predictions = ['Unknown'] * n_cells
    confidence_scores = [0.0] * n_cells
    
    return predictions, confidence_scores

def main():
    # Set paths
    SEURAT_PATH = "../data/Liu_Dataset.RDS"
    output_dir = "scMatch_output"
    
    print("=== scMatch Benchmarking ===")
    
    start_total = time.time()
    
    try:
        # Load data
        expr_df, metadata_df = load_seurat_object(SEURAT_PATH)
        
        print(f"Loaded data: {expr_df.shape[1]} cells x {expr_df.shape[0]} genes")
        print(f"Cell types: {metadata_df['Ground_Truth_Celltype'].nunique()} unique types")
        
        # Perform 70/30 split
        train_expr, test_expr, train_meta, test_meta = stratified_split(expr_df, metadata_df)
        
        # Create custom reference from training data
        reference_dir = create_custom_reference(train_expr, train_meta, output_dir)
        
        # Run scMatch
        start_time = time.time()
        
        scmatch_output = run_scmatch(test_expr, reference_dir, output_dir)
        
        if scmatch_output and os.path.exists(scmatch_output):
            # Parse results
            predictions, cell_ids = parse_scmatch_results(scmatch_output, test_meta)
            
            if predictions is not None:
                # Create confidence scores (scMatch may not provide these directly)
                confidence_scores = [1.0 if pred != 'Unknown' else 0.0 for pred in predictions]
            else:
                predictions, confidence_scores = create_dummy_results(test_meta)
        else:
            predictions, confidence_scores = create_dummy_results(test_meta)
        
        end_time = time.time()
        runtime = end_time - start_time
        
        # Create standardized output
        output_df = pd.DataFrame({
            'cell_id': test_meta.index,
            'predicted_type': predictions,
            'confidence_score': confidence_scores,
            'true_type': test_meta['Ground_Truth_Celltype'].values
        })
        
        # Calculate metrics
        if not all(pred == 'Unknown' for pred in predictions):
            metrics = calculate_metrics(output_df['predicted_type'], output_df['true_type'])
            
            # Display results
            total_cells = len(output_df)
            correct_cells = sum(output_df['predicted_type'] == output_df['true_type'])
            
            print(f"\\n=== scMatch Results Summary ===")
            print(f"Runtime: {runtime:.2f} seconds")
            print(f"Total cells processed: {total_cells}\\n")
            
            print("Performance Metrics:")
            print(f"- Overall Accuracy: {metrics['overall_accuracy']*100:.2f}% ({correct_cells}/{total_cells} cells correct)")
            print(f"- Macro F1 Score: {metrics['macro_f1']:.3f}\\n")
            
            print("Per-Class Performance:")
            per_class_df = metrics['per_class_metrics']
            per_class_df = per_class_df.round({'precision': 3, 'recall': 3, 'f1': 3})
            print(per_class_df.to_string(index=False))
            
            # Save metrics summary
            per_class_df.to_csv("scMatch_metrics_summary.csv", index=False)
            print("\\nPer-class metrics saved to scMatch_metrics_summary.csv")
            
            # Confusion matrix
            from sklearn.metrics import confusion_matrix
            unique_labels = np.union1d(output_df['true_type'].unique(), output_df['predicted_type'].unique())
            cm = confusion_matrix(output_df['true_type'], output_df['predicted_type'], labels=unique_labels)
            
            print("\\nConfusion Matrix:")
            cm_df = pd.DataFrame(cm, index=unique_labels, columns=unique_labels)
            print(cm_df)
        else:
            print("\\n=== scMatch Results Summary ===")
            print(f"Runtime: {runtime:.2f} seconds")
            print(f"Total cells processed: {len(output_df)}")
            print("All predictions were 'Unknown' - scMatch may have failed")
        
        print("\\nPredicted cell types distribution:")
        print(output_df['predicted_type'].value_counts())
        
        # Save results
        output_df.to_csv("scMatch_results.csv", index=False)
        print("\\nResults saved to scMatch_results.csv")
        
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        print("Creating minimal output file...")
        
        # Create minimal dummy output
        dummy_output = pd.DataFrame({
            'cell_id': ['cell_1'],
            'predicted_type': ['Unknown'],
            'confidence_score': [0.0],
            'true_type': ['Unknown']
        })
        dummy_output.to_csv("scMatch_results.csv", index=False)
        print("Dummy results saved to scMatch_results.csv")
    
    total_runtime = time.time() - start_total
    print(f"\\nTotal execution time: {total_runtime:.2f} seconds")

if __name__ == "__main__":
    main()