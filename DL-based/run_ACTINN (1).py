# run_ACTINN.py
#################################################
# ACTINN Function for Python Benchmarking Framework
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
import h5py
from typing import Dict, List, Any


def run_ACTINN_function(adata_train: ad.AnnData, adata_test: ad.AnnData, markers: pd.DataFrame) -> Dict[str, Any]:
    """
    ACTINN Cell Type Annotation Function

    Purpose: Run ACTINN algorithm using a 4-layer neural network for cell type classification
    Inputs:
      - adata_train: Training AnnData object (used for model training)
      - adata_test: Test AnnData object to predict
      - markers: Marker genes dataframe from rank_genes_groups() (for interface consistency)
    Outputs: Dict with predictions, true_labels, confidence_scores, cell_ids
    Algorithm: Uses ACTINN's 4-layer fully connected neural network (100->50->25->n_types)
               with Adam optimizer and L2 regularization
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
    conda_env_path = "/home/oliver/miniconda3/envs/ACTINN_env"
    if not os.path.isdir(conda_env_path):
        warnings.warn(f"ACTINN_env conda environment not found at {conda_env_path}. Please create it first.")
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

            # Save AnnData objects for inter-process communication
            # Strip to minimal data for old anndata (0.7.6) compatibility
            for adata in [adata_train_subset, adata_test_subset]:
                adata.obs = adata.obs[['Ground_Truth_Celltype']].copy()
                adata.obs['Ground_Truth_Celltype'] = adata.obs['Ground_Truth_Celltype'].astype(str)
                adata.var = adata.var[[]]
                adata.uns.clear()
                adata.layers.clear()
            adata_train_subset.write_h5ad(train_h5ad_path)
            adata_test_subset.write_h5ad(test_h5ad_path)
            # Patch h5ad files for old anndata (0.7.6) compatibility in ACTINN_env
            for h5_path in [train_h5ad_path, test_h5ad_path]:
                with h5py.File(h5_path, 'a') as f:
                    # Remove /layers group — old anndata can't parse encoding-type 'dict'
                    if 'layers' in f:
                        del f['layers']
                    # Convert categorical columns in /obs to plain string datasets
                    # New anndata writes categoricals as groups with 'codes' + 'categories'
                    # but old anndata/h5py can't read them
                    if 'obs' in f:
                        for col_name in list(f['obs'].keys()):
                            col = f['obs'][col_name]
                            if isinstance(col, h5py.Group) and 'codes' in col and 'categories' in col:
                                codes = col['codes'][()]
                                categories = col['categories'][()]
                                if hasattr(categories[0], 'decode'):
                                    categories = [c.decode('utf-8') for c in categories]
                                string_values = [categories[c] for c in codes]
                                del f['obs'][col_name]
                                f['obs'].create_dataset(col_name, data=[s.encode('utf-8') for s in string_values])

            # Create Python script for ACTINN execution with embedded algorithm
            actinn_script = f'''
import tracemalloc
tracemalloc.start()  # Start tracking memory

import sys
import os
import numpy as np
import pandas as pd
import math
import pickle
import warnings
warnings.filterwarnings("ignore")

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
from tensorflow.python.framework import ops
import scipy.sparse
import scanpy as sc

# Set seeds for reproducibility
np.random.seed(1)
tf.set_random_seed(3)

# ============================================
# ACTINN Core Functions (embedded from actinn_predict.py)
# ============================================

def adata_to_actinn_df(adata):
    """Convert AnnData to DataFrame format ACTINN expects (genes x cells, uppercase)"""
    if scipy.sparse.issparse(adata.X):
        data = adata.X.T.toarray()  # genes x cells
    else:
        data = adata.X.T
    return pd.DataFrame(
        data,
        index=[g.upper() for g in adata.var_names],
        columns=adata.obs_names
    )

def scale_sets(sets):
    """Get common genes, normalize and scale the sets"""
    common_genes = set(sets[0].index)
    for i in range(1, len(sets)):
        common_genes = set.intersection(set(sets[i].index), common_genes)
    common_genes = sorted(list(common_genes))
    sep_point = [0]
    for i in range(len(sets)):
        sets[i] = sets[i].loc[common_genes,]
        sep_point.append(sets[i].shape[1])
    total_set = np.array(pd.concat(sets, axis=1, sort=False), dtype=np.float32)
    total_set = np.divide(total_set, np.sum(total_set, axis=0, keepdims=True)) * 10000
    total_set = np.log2(total_set + 1)
    expr = np.sum(total_set, axis=1)
    total_set = total_set[np.logical_and(expr >= np.percentile(expr, 1), expr <= np.percentile(expr, 99)),]
    cv = np.std(total_set, axis=1) / np.mean(total_set, axis=1)
    total_set = total_set[np.logical_and(cv >= np.percentile(cv, 1), cv <= np.percentile(cv, 99)),]
    for i in range(len(sets)):
        sets[i] = total_set[:, sum(sep_point[:(i+1)]):sum(sep_point[:(i+2)])]
    return sets

def one_hot_matrix(labels, C):
    """Turn labels into one-hot matrix"""
    C = tf.constant(C, name="C")
    one_hot_matrix = tf.one_hot(labels, C, axis=0)
    sess = tf.Session()
    one_hot = sess.run(one_hot_matrix)
    sess.close()
    return one_hot

def type_to_label_dict(types):
    """Make types to labels dictionary"""
    type_to_label = {{}}
    all_type = list(set(types))
    for i in range(len(all_type)):
        type_to_label[all_type[i]] = i
    return type_to_label

def convert_type_to_label(types, type_to_label):
    """Convert types to labels"""
    types = list(types)
    labels = []
    for t in types:
        labels.append(type_to_label[t])
    return labels

def create_placeholders(n_x, n_y):
    """Create TensorFlow placeholders"""
    X = tf.placeholder(tf.float32, shape=(n_x, None))
    Y = tf.placeholder(tf.float32, shape=(n_y, None))
    return X, Y

def initialize_parameters(nf, ln1, ln2, ln3, nt):
    """Initialize parameters with Xavier initialization"""
    tf.set_random_seed(3)
    W1 = tf.get_variable("W1", [ln1, nf], initializer=tf.contrib.layers.xavier_initializer(seed=3))
    b1 = tf.get_variable("b1", [ln1, 1], initializer=tf.zeros_initializer())
    W2 = tf.get_variable("W2", [ln2, ln1], initializer=tf.contrib.layers.xavier_initializer(seed=3))
    b2 = tf.get_variable("b2", [ln2, 1], initializer=tf.zeros_initializer())
    W3 = tf.get_variable("W3", [ln3, ln2], initializer=tf.contrib.layers.xavier_initializer(seed=3))
    b3 = tf.get_variable("b3", [ln3, 1], initializer=tf.zeros_initializer())
    W4 = tf.get_variable("W4", [nt, ln3], initializer=tf.contrib.layers.xavier_initializer(seed=3))
    b4 = tf.get_variable("b4", [nt, 1], initializer=tf.zeros_initializer())
    parameters = {{"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}}
    return parameters

def forward_propagation(X, parameters):
    """Forward propagation: LINEAR->RELU->LINEAR->RELU->LINEAR->RELU->LINEAR"""
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4

def forward_propagation_for_predict(X, parameters):
    """Forward propagation for prediction"""
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4

def compute_cost(Z4, Y, parameters, lambd=0.01):
    """Compute cost with L2 regularization"""
    logits = tf.transpose(Z4)
    labels = tf.transpose(Y)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=logits, labels=labels)) + \\
        (tf.nn.l2_loss(parameters["W1"]) + tf.nn.l2_loss(parameters["W2"]) +
         tf.nn.l2_loss(parameters["W3"]) + tf.nn.l2_loss(parameters["W4"])) * lambd
    return cost

def random_mini_batches(X, Y, mini_batch_size=32, seed=1):
    """Generate random mini batches"""
    ns = X.shape[1]
    mini_batches = []
    np.random.seed(seed)
    permutation = list(np.random.permutation(ns))
    shuffled_X = X[:, permutation]
    shuffled_Y = Y[:, permutation]
    num_complete_minibatches = int(math.floor(ns / mini_batch_size))
    for k in range(0, num_complete_minibatches):
        mini_batch_X = shuffled_X[:, k * mini_batch_size: k * mini_batch_size + mini_batch_size]
        mini_batch_Y = shuffled_Y[:, k * mini_batch_size: k * mini_batch_size + mini_batch_size]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    if ns % mini_batch_size != 0:
        mini_batch_X = shuffled_X[:, num_complete_minibatches * mini_batch_size: ns]
        mini_batch_Y = shuffled_Y[:, num_complete_minibatches * mini_batch_size: ns]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    return mini_batches

def predict(X, parameters):
    """Get predicted class indices"""
    W1 = tf.convert_to_tensor(parameters["W1"])
    b1 = tf.convert_to_tensor(parameters["b1"])
    W2 = tf.convert_to_tensor(parameters["W2"])
    b2 = tf.convert_to_tensor(parameters["b2"])
    W3 = tf.convert_to_tensor(parameters["W3"])
    b3 = tf.convert_to_tensor(parameters["b3"])
    W4 = tf.convert_to_tensor(parameters["W4"])
    b4 = tf.convert_to_tensor(parameters["b4"])
    params = {{"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}}
    x = tf.placeholder("float")
    z4 = forward_propagation_for_predict(x, params)
    p = tf.argmax(z4)
    sess = tf.Session()
    prediction = sess.run(p, feed_dict={{x: X}})
    sess.close()
    return prediction

def predict_probability(X, parameters):
    """Get softmax probabilities for each class"""
    W1 = tf.convert_to_tensor(parameters["W1"])
    b1 = tf.convert_to_tensor(parameters["b1"])
    W2 = tf.convert_to_tensor(parameters["W2"])
    b2 = tf.convert_to_tensor(parameters["b2"])
    W3 = tf.convert_to_tensor(parameters["W3"])
    b3 = tf.convert_to_tensor(parameters["b3"])
    W4 = tf.convert_to_tensor(parameters["W4"])
    b4 = tf.convert_to_tensor(parameters["b4"])
    params = {{"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}}
    x = tf.placeholder("float")
    z4 = forward_propagation_for_predict(x, params)
    p = tf.nn.softmax(z4, axis=0)
    sess = tf.Session()
    prediction = sess.run(p, feed_dict={{x: X}})
    sess.close()
    return prediction

def model(X_train, Y_train, X_test, starting_learning_rate=0.0001, num_epochs=50, minibatch_size=128, print_cost=True):
    """Build and train the 4-layer neural network model"""
    ops.reset_default_graph()
    tf.set_random_seed(3)
    seed = 3
    (nf, ns) = X_train.shape
    nt = Y_train.shape[0]
    costs = []
    X, Y = create_placeholders(nf, nt)
    parameters = initialize_parameters(nf=nf, ln1=100, ln2=50, ln3=25, nt=nt)
    Z4 = forward_propagation(X, parameters)
    cost = compute_cost(Z4, Y, parameters, 0.005)
    global_step = tf.Variable(0, trainable=False)
    learning_rate = tf.train.exponential_decay(starting_learning_rate, global_step, 1000, 0.95, staircase=True)
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate)
    trainer = optimizer.minimize(cost, global_step=global_step)
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        for epoch in range(num_epochs):
            epoch_cost = 0.
            num_minibatches = int(ns / minibatch_size)
            seed = seed + 1
            minibatches = random_mini_batches(X_train, Y_train, minibatch_size, seed)
            for minibatch in minibatches:
                (minibatch_X, minibatch_Y) = minibatch
                _, minibatch_cost = sess.run([trainer, cost], feed_dict={{X: minibatch_X, Y: minibatch_Y}})
                epoch_cost += minibatch_cost / num_minibatches
            if print_cost and (epoch + 1) % 5 == 0:
                costs.append(epoch_cost)
        parameters = sess.run(parameters)
        return parameters

# ============================================
# Main Execution
# ============================================
try:
    adata_train = sc.read_h5ad('{train_h5ad_path}')
    adata_test = sc.read_h5ad('{test_h5ad_path}')

    # Convert to ACTINN DataFrame format (genes x cells, uppercase)
    train_set = adata_to_actinn_df(adata_train)
    test_set = adata_to_actinn_df(adata_test)

    # Remove duplicate genes (keep first occurrence)
    train_set = train_set.loc[~train_set.index.duplicated(keep='first')]
    test_set = test_set.loc[~test_set.index.duplicated(keep='first')]

    # Get labels from training data
    train_labels_raw = adata_train.obs['Ground_Truth_Celltype'].tolist()
    barcode = list(test_set.columns)
    nt = len(set(train_labels_raw))

    # Scale datasets
    train_set, test_set = scale_sets([train_set, test_set])

    # Convert labels to numeric
    type_to_label = type_to_label_dict(train_labels_raw)
    label_to_type = {{v: k for k, v in type_to_label.items()}}
    train_labels_numeric = convert_type_to_label(train_labels_raw, type_to_label)
    train_labels_onehot = one_hot_matrix(train_labels_numeric, nt)

    # Train model
    parameters = model(
        train_set,
        train_labels_onehot,
        test_set,
        starting_learning_rate=0.0001,
        num_epochs=50,
        minibatch_size=128,
        print_cost=True
    )

    # Get predictions
    test_predict_indices = predict(test_set, parameters)
    predictions = [label_to_type[idx] for idx in test_predict_indices]

    # Get confidence scores from probability matrix
    prob_matrix = predict_probability(test_set, parameters)
    confidence_scores = np.max(prob_matrix, axis=0).tolist()

    # Get true labels
    true_labels = adata_test.obs['Ground_Truth_Celltype'].tolist()
    cell_ids = list(adata_test.obs.index)

    # Prepare results
    results = {{
        'predictions': predictions,
        'confidence_scores': confidence_scores,
        'true_labels': true_labels,
        'cell_ids': cell_ids
    }}

    # Capture peak memory usage
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results['peak_memory_mb'] = peak_mem / (1024 * 1024)

    print(f"Peak memory usage: {{results['peak_memory_mb']:.2f}} MB")

    # Save results
    with open('{results_path}', 'wb') as f:
        pickle.dump(results, f)

    print(f"Predictions: {{len(predictions)}} cells, unique: {{set(predictions)}}")

except Exception as e:
    print(f"ACTINN execution failed: {{str(e)}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''

            # Execute ACTINN in conda environment
            # Set environment - enable GPU
            env = os.environ.copy()
            env['MPLBACKEND'] = 'Agg'  # Use non-interactive backend

            result = subprocess.run(
                [
                    "/home/oliver/miniconda3/condabin/conda",
                    "run",
                    "-n", "ACTINN_env",
                    "python",
                    "-c", actinn_script
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=86400  # 24 hour timeout for training
            )

            if result.returncode != 0:
                warnings.warn(f"ACTINN execution failed: {result.stderr}")
                print(f"stdout: {result.stdout}")
                return default_return()

            # Print subprocess output (contains verbose prediction comparison)
            print(result.stdout)

            # Load results from pickle file
            if not os.path.exists(results_path):
                warnings.warn("ACTINN results file not found")
                return default_return()

            with open(results_path, 'rb') as f:
                results = pickle.load(f)

            cell_predictions = results['predictions']
            confidence_scores = results['confidence_scores']
            true_labels = results['true_labels']
            cell_ids = results['cell_ids']
            peak_memory_mb = results.get('peak_memory_mb', None)

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

            print(f"ACTINN completed successfully:")
            print(f"  - Total predictions: {len(cell_predictions)}")
            print(f"  - Unknown predictions: {sum(1 for p in cell_predictions if p == 'Unknown')}")
            print(f"  - Unique predictions: {len(set(cell_predictions))}")

            # Return standardized format
            return {
                'predictions': cell_predictions,
                'true_labels': true_labels,
                'confidence_scores': confidence_scores,
                'cell_ids': cell_ids,
                'peak_memory_mb': peak_memory_mb
            }

        except subprocess.TimeoutExpired:
            warnings.warn("ACTINN execution timed out after 24 hours")
            return default_return()
        except Exception as e:
            warnings.warn(f"ACTINN error: {str(e)}")
            return default_return()


# For backward compatibility
run_ACTINN = run_ACTINN_function
