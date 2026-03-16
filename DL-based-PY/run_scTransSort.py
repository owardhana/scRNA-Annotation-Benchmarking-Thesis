# run_scTransSort.py
#################################################
# scTransSort Function for Python Benchmarking Framework
# Input: Train/test AnnData objects and markers from rank_genes_groups
# Output: Standardized results format for CV framework
#################################################

import subprocess
import tempfile
import os
import shutil
import pandas as pd
import anndata as ad
import warnings
import pickle
from typing import Dict, List, Any


def run_scTransSort_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    scTransSort Cell Type Annotation Function

    Purpose: Run scTransSort algorithm using Vision Transformer for cell type classification
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses scTransSort's Vision Transformer with gene expression as image embeddings

    Reference:
    Jiao et al. (2023). scTransSort: Transformers for Intelligent Annotation of Cell Types
    by Gene Embeddings. Biomolecules, 13(4), 611.
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
    conda_env_path = "/home/oliver/miniconda3/envs/scTransSort_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"scTransSort_env conda environment not found at {conda_env_path}. Please create it first.")
        return default_return()

    # Create temporary directory for file operations
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Preparing data for scTransSort...")

            # Prepare file paths
            train_h5ad_path = os.path.join(temp_dir, "adata_train.h5ad")
            test_h5ad_path = os.path.join(temp_dir, "adata_test.h5ad")
            model_weights_path = os.path.join(temp_dir, "scTransSort_model.weights.h5")
            results_path = os.path.join(temp_dir, "results.pkl")

            # Get common genes between train and test
            common_genes = list(set(adata_train.var.index) & set(adata_test.var.index))
            if len(common_genes) < 100:
                warnings.warn("Too few common genes between train and test datasets")
                return default_return()

            print(f"Using {len(common_genes)} common genes for scTransSort analysis")

            # Subset both datasets to common genes and ensure consistent ordering
            adata_train_subset = adata_train[:, sorted(common_genes)].copy()
            adata_test_subset = adata_test[:, sorted(common_genes)].copy()

            # Save AnnData objects for inter-process communication
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)

            # Copy scTransSort_core module to temp directory for subprocess access
            core_src = os.path.join(os.path.dirname(__file__), "scTransSort_core")
            core_dst = os.path.join(temp_dir, "scTransSort_core")
            shutil.copytree(core_src, core_dst)

            print("Running scTransSort in conda environment...")

            # Create Python script for scTransSort execution
            sctranssort_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory (before imports, after data conversion)

import sys
import os
import pickle
import warnings
warnings.filterwarnings("ignore")

# CRITICAL: Configure TensorFlow BEFORE importing it
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Suppress TF warnings (3 = ERROR only)
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'  # Enable memory growth from start

# Force CPU mode for problematic GPUs (RTX 50xx series with compute capability 12.0+)
# These GPUs are too new for most TensorFlow builds and cause JIT compilation issues
FORCE_CPU = True  # Set to False to attempt GPU (not recommended for RTX 5060)
USE_GPU = False  # Will be set based on GPU detection
GPU_DEVICE = '/GPU:0'
CPU_DEVICE = '/CPU:0'

try:
    import tensorflow as tf

    if FORCE_CPU:
        print("⚠ CPU-only mode enabled (FORCE_CPU=True)")
        print("  Reason: RTX 5060 compute capability 12.0 incompatible with TensorFlow")
        print("  Training will use CPU (slower but stable)")
        USE_GPU = False
    else:
        # Attempt GPU detection and testing
        gpus = tf.config.list_physical_devices('GPU')
        if gpus:
            try:
                # Enable memory growth
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)

                # Check compute capability
                gpu_details = tf.config.experimental.get_device_details(gpus[0])
                compute_cap = gpu_details.get('compute_capability', (0, 0))
                compute_version = compute_cap[0] + compute_cap[1] / 10.0

                print(f"✓ Found {{len(gpus)}} GPU(s)")
                print(f"  GPU: {{gpu_details.get('device_name', 'Unknown')}}")
                print(f"  Compute capability: {{compute_cap[0]}}.{{compute_cap[1]}}")

                # Check if compute capability is too new (>= 12.0)
                if compute_version >= 12.0:
                    print(f"⚠ Compute capability {{compute_version}} too new for TensorFlow build")
                    print(f"⚠ Would require JIT compilation (30+ min), forcing CPU instead")
                    USE_GPU = False
                else:
                    # Test GPU with a simple operation
                    with tf.device(GPU_DEVICE):
                        test_tensor = tf.constant([1.0, 2.0, 3.0])
                        test_result = tf.reduce_sum(test_tensor)
                    print(f"✓ GPU test passed - will use GPU for training")
                    USE_GPU = True

            except Exception as e:
                print(f"⚠ GPU initialization failed: {{str(e)[:100]}}")
                print(f"⚠ Falling back to CPU")
                USE_GPU = False
        else:
            print("⚠ No GPU detected - training will use CPU")
            USE_GPU = False

    import numpy as np
    import anndata as ad
    from collections import Counter

    # Import scTransSort core modules
    sys.path.insert(0, '{temp_dir}')
    from scTransSort_core.data_processor import (
        reshape_to_image, prepare_tf_dataset, encode_celltypes, clean_expression_data
    )
    from scTransSort_core.vit_model import build_sctranssort_model, create_cosine_lr_schedule
    from scTransSort_core.config import ScTransSortConfig

    print("=" * 60)
    print("scTransSort Training and Prediction")
    print("=" * 60)

    # Load data
    print("\\nLoading training and test data...")
    adata_train = ad.read_h5ad('{train_h5ad_path}')
    adata_test = ad.read_h5ad('{test_h5ad_path}')

    print(f"Training data: {{len(adata_train.obs)}} cells x {{len(adata_train.var)}} genes")
    print(f"Test data: {{len(adata_test.obs)}} cells x {{len(adata_test.var)}} genes")

    # Extract and encode cell types
    train_celltypes = adata_train.obs['Ground_Truth_Celltype'].values
    test_celltypes = adata_test.obs['Ground_Truth_Celltype'].values

    train_labels_int, num_classes, unique_celltypes = encode_celltypes(train_celltypes)

    print(f"\\nCell type statistics:")
    print(f"  Number of unique cell types: {{num_classes}}")
    print(f"  Cell types: {{list(unique_celltypes)}}")

    train_dist = Counter(train_celltypes)
    print(f"\\nTraining cell type distribution:")
    for celltype, count in sorted(train_dist.items()):
        print(f"  {{celltype}}: {{count}} cells")

    # Clean expression data
    print("\\nCleaning expression data...")
    train_expr_clean = clean_expression_data(adata_train.X, verbose=True)
    test_expr_clean = clean_expression_data(adata_test.X, verbose=True)

    # Reshape to image format (224x224x3)
    print("\\nReshaping data to image format...")
    train_images = reshape_to_image(train_expr_clean, target_size=224)
    test_images = reshape_to_image(test_expr_clean, target_size=224)

    # Create configuration
    config = ScTransSortConfig(
        img_size=224,
        patch_size=16,
        embed_dim=768,
        depth=12,
        num_heads=12,
        batch_size=64,
        epochs=40,
        learning_rate=0.001,
        val_split=0.3
    )

    print(f"\\n{{config}}")

    # Create TensorFlow datasets
    print("\\nCreating optimized TensorFlow datasets...")
    train_ds, val_ds = prepare_tf_dataset(
        train_images,
        train_labels_int,
        batch_size=config.batch_size,
        val_split=config.val_split,
        shuffle=True,
        cache_data=False,  # Set to True if sufficient RAM
        img_size=config.img_size
    )

    # Build model (with device context and error handling)
    print("\\nBuilding scTransSort model...")
    device = GPU_DEVICE if USE_GPU else CPU_DEVICE
    print(f"Using device: {{device}}")

    try:
        with tf.device(device):
            model = build_sctranssort_model(num_classes=num_classes, config=config)
    except Exception as e:
        if USE_GPU:
            print(f"\\n⚠ Model building failed on GPU: {{str(e)[:100]}}")
            print(f"⚠ Retrying with CPU...")
            USE_GPU = False
            device = CPU_DEVICE
            with tf.device(device):
                model = build_sctranssort_model(num_classes=num_classes, config=config)
        else:
            raise  # Re-raise if already on CPU

    # Compile model
    print("\\nCompiling model...")

    # Loss function
    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

    # Optimizer
    optimizer = tf.keras.optimizers.SGD(
        learning_rate=config.learning_rate,
        momentum=config.momentum
    )

    # Metrics
    metrics = [
        tf.keras.metrics.SparseCategoricalAccuracy(name='accuracy')
    ]

    model.compile(
        optimizer=optimizer,
        loss=loss_fn,
        metrics=metrics
    )

    # Learning rate scheduler
    lr_scheduler = tf.keras.callbacks.LearningRateScheduler(
        create_cosine_lr_schedule(config.learning_rate, config.epochs, config.end_lr_rate)
    )

    # Model checkpoint (save best model)
    checkpoint_cb = tf.keras.callbacks.ModelCheckpoint(
        '{model_weights_path}',
        monitor='val_accuracy',
        save_best_only=True,
        save_weights_only=True,
        mode='max',
        verbose=1
    )

    # Early stopping (optional, disabled by default)
    # early_stop_cb = tf.keras.callbacks.EarlyStopping(
    #     monitor='val_accuracy',
    #     patience=10,
    #     restore_best_weights=True
    # )

    callbacks = [lr_scheduler, checkpoint_cb]

    # Train model
    print("\\n" + "=" * 60)
    print("Training scTransSort model...")
    print("=" * 60)

    history = model.fit(
        train_ds,
        validation_data=val_ds,
        epochs=config.epochs,
        callbacks=callbacks,
        verbose=1
    )

    print("\\nTraining completed!")

    # Load best model weights
    print("\\nLoading best model weights...")
    model.load_weights('{model_weights_path}')

    # Make predictions on test set
    print("\\nMaking predictions on test set...")

    # Encode test labels for comparison
    test_celltype_to_int = {{ct: i for i, ct in enumerate(unique_celltypes)}}
    test_labels_int = np.array([
        test_celltype_to_int.get(ct, -1) for ct in test_celltypes
    ])

    # Create test dataset (no shuffling, no splitting)
    test_ds = tf.data.Dataset.from_tensor_slices((test_images, test_labels_int))
    test_ds = test_ds.batch(config.batch_size).prefetch(tf.data.AUTOTUNE)

    # Predict
    predictions_logits = model.predict(test_ds, verbose=1)

    # Convert logits to probabilities
    predictions_probs = tf.nn.softmax(predictions_logits, axis=-1).numpy()

    # Get predicted class indices
    predicted_indices = np.argmax(predictions_probs, axis=-1)

    # Get confidence scores (max probability for each prediction)
    confidence_scores = np.max(predictions_probs, axis=-1)

    # Map indices back to cell type names
    predicted_celltypes = [unique_celltypes[idx] for idx in predicted_indices]
    true_celltypes = list(test_celltypes)

    print(f"\\nPrediction completed:")
    print(f"  Total predictions: {{len(predicted_celltypes)}}")
    print(f"  Unique predicted types: {{len(set(predicted_celltypes))}}")

    # Calculate accuracy
    correct = sum([p == t for p, t in zip(predicted_celltypes, true_celltypes)])
    accuracy = correct / len(true_celltypes)
    print(f"  Test accuracy: {{accuracy:.4f}} ({{correct}}/{{len(true_celltypes)}})")

    # Prediction statistics
    print(f"\\n=== Prediction vs Ground Truth Comparison ===")
    pred_counts = Counter(predicted_celltypes)
    true_counts = Counter(true_celltypes)

    print(f"\\nGround Truth Distribution:")
    for celltype, count in sorted(true_counts.items()):
        print(f"  {{celltype}}: {{count}} cells")

    print(f"\\nPredicted Distribution:")
    for celltype, count in sorted(pred_counts.items()):
        print(f"  {{celltype}}: {{count}} cells")

    # Per-class metrics
    print(f"\\nPer-Class Accuracy:")
    from sklearn.metrics import classification_report
    print(classification_report(true_celltypes, predicted_celltypes, zero_division=0))

    # Show sample predictions
    cell_ids = list(adata_test.obs.index)
    print(f"\\nSample Predictions (first 20 cells):")
    print(f"{{'Cell ID':<20}} {{'Predicted':<20}} {{'Actual':<20}} {{'Confidence':<12}} {{'Correct':<10}}")
    print("-" * 82)
    for i in range(min(20, len(predicted_celltypes))):
        correct_mark = "✓" if predicted_celltypes[i] == true_celltypes[i] else "✗"
        print(f"{{cell_ids[i]:<20}} {{predicted_celltypes[i]:<20}} {{true_celltypes[i]:<20}} {{confidence_scores[i]:<12.4f}} {{correct_mark:<10}}")

    # Prepare results
    results = {{
        'predictions': predicted_celltypes,
        'confidence_scores': list(confidence_scores),
        'true_labels': true_celltypes,
        'cell_ids': cell_ids,
        'accuracy': accuracy,
        'num_classes': num_classes,
        'unique_celltypes': list(unique_celltypes)
    }}

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)  # Convert bytes to MB

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print("\\n" + "=" * 60)
    print("scTransSort execution completed successfully!")
    print("=" * 60)

except Exception as e:
    print(f"scTransSort execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute scTransSort in conda environment
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend
            env['PYTHONPATH'] = temp_dir  # Add temp dir to Python path

            # Force CPU mode by hiding GPU from TensorFlow
            # This prevents the RTX 5060 compute capability 12.0 incompatibility
            env['CUDA_VISIBLE_DEVICES'] = ''  # Hide all GPUs from TensorFlow
            print("Note: Running in CPU-only mode (CUDA_VISIBLE_DEVICES='')")
            print("      Reason: RTX 5060 not compatible with current TensorFlow build")

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "scTransSort_env",
                    "python",
                    "-c", sctranssort_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour (1 day) timeout
            )

            if result.returncode != 0:
                warnings.warn(f"scTransSort execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains training progress and results)
            print(result.stdout)

            print("Parsing scTransSort results...")

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("scTransSort results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidence_scores = results['confidence_scores']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            test_accuracy = results.get('accuracy', 0.0)
            peak_memory_mb = results.get('peak_memory_mb', None)  # Extract peak memory

            # Ensure predictions have correct length
            if len(cell_predictions) != len(adata_test.obs):
                warnings.warn(f"Prediction length mismatch. Expected: {len(adata_test.obs)}, Got: {len(cell_predictions)}")
                return default_return()

            # Convert predictions to strings and handle missing values
            cell_predictions = [str(pred) if pd.notna(pred) else "Unknown" for pred in cell_predictions]

            # Ensure confidence scores are valid floats in [0, 1]
            confidence_scores = [
                max(0.0, min(1.0, float(conf))) if pd.notna(conf) and isinstance(conf, (int, float)) else 0.5
                for conf in confidence_scores
            ]

            print(f"scTransSort completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")
            print(f"  - Test accuracy: {test_accuracy:.4f}")
            print(f"  - Mean confidence: {np.mean(confidence_scores):.4f}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb  # Add peak memory to return dict
            }

        except subprocess.TimeoutExpired:
            warnings.warn("scTransSort execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"scTransSort error: {str(e)}")
            import traceback
            traceback.print_exc()
            return default_return()


# For backward compatibility
run_scTransSort = run_scTransSort_function
