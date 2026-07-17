# callr_core.R
#################################################
# CALLR Core Functions
# Source: https://github.com/MathSZhang/CALLR
# Cell type Annotation using Laplacian and Logistic Regression
#################################################

# Load required libraries
library(glmnet)
library(Matrix)

#' Distance computation function
#' Computes distance matrix between data points
dist2 <- function(x, c = NA) {

  # set the parameters for x
  if(is.na(c)) {
    c = x
  }

  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }

  # compute the distance
  dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) +
    (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) -
    2 * (x%*%t(c))

  return(dist)
}

#' Simplex projection function
#' Projects a vector to a probability simplex
ps <- function(y){

  z=sort(y)
  n=length(y)
  t=seq(n)
  for(i in 1:n-1){
    t[i]=(sum(z[(i+1):n])-1)/(n-i)
  }
  t[n]=z[n]-1
  tt=n-length(which(t<z))

  if(tt==0){
    tt=(sum(y)-1)/n
  }else{
    tt=t[tt]
  }

  x=y-tt
  x[x<0]<-0
  return(x)
}

#' Training set generation function
#' Creates a training set matrix from a label vector
rand <- function(x, r){

  aaa=c()
  for (i in 1:length(unique(x))){
    aaa=c(aaa,sample(which(x==sort(unique(x))[i]),
                     max(ceiling(length(which(x==sort(unique(x))[i]))*r),2),replace=FALSE))
  }
  F=matrix(0, length(x), length(unique(x)))
  for(i in 1:length(unique(x))){
    F[which(x==sort(unique(x))[i]),i] = 1
  }

  return(list("training_set" = cbind(aaa, x[aaa]), "F" = F))
}

#' Main CALLR classification function
#'
#' Implements the original MBO-based semi-supervised algorithm from Zhang et al.
#' Two nested loops:
#'   Outer: binarize U -> refit LR on all cells -> update U_hat -> re-pin training labels
#'   Inner: MBO gradient steps (-L*U + u*U_hat) with simplex projection per row
#'
#' @param X  Numeric matrix, cells x features (samples x genes)
#' @param u  Scalar, weight for the logistic regression signal (original default ~0.5)
#' @param T_in  2-column matrix: col1 = cell row indices (1-based), col2 = integer class labels
#' @return Integer vector of length nrow(X), predicted class index for each cell (1-based, aligned with sort(unique(T_in[,2])))
callr <- function(X, u, T_in) {

  n  <- nrow(X)
  nc_raw <- length(unique(T_in[,2]))

  # --- Normalise class labels to contiguous 1:nc to avoid factor-level gaps ---
  # (factor("10") sorts before factor("2") alphabetically, breaking glmnet column order)
  label_vals <- sort(unique(T_in[,2]))  # sorted numeric class IDs
  T <- cbind(T_in[,1], match(T_in[,2], label_vals))  # remap labels to 1:nc
  nc <- length(label_vals)

  # --- Compute distance + weight matrix ---
  D <- dist2(X)

  k_neighbors <- min(10, max(2, floor(n / 4)))
  sigma <- apply(D, 1, function(row) {
    s <- sort(row[row > 0])
    if (length(s) >= k_neighbors) s[k_neighbors]
    else if (length(s) > 0)       s[length(s)]
    else                           1
  })
  sigma <- pmax(sigma, quantile(sigma, 0.1, na.rm = TRUE))

  W <- exp(-D / outer(sigma, sigma))
  diag(W) <- 0

  # --- Normalised symmetric Laplacian: L = I - D^{-1/2} W D^{-1/2} ---
  # (original uses normalised Laplacian, not combinatorial D - W)
  d <- rowSums(W)
  d[d == 0] <- 1                   # guard against isolated nodes
  d_inv_sqrt <- 1 / sqrt(d)
  L <- diag(n) - W * outer(d_inv_sqrt, d_inv_sqrt)

  # --- Helper: pin training rows of U to one-hot ---
  pin_training <- function(U_mat) {
    for (k in 1:nrow(T)) {
      U_mat[T[k,1], ]     <- 0
      U_mat[T[k,1], T[k,2]] <- 1
    }
    return(U_mat)
  }

  # --- Helper: fit LR and return U_hat (log-prob matrix, training rows hard-pinned) ---
  #
  # x_fit: features for training observations
  # y_fit: integer labels for training observations (in 1:nc)
  # X_all: features for all n cells (used for prediction)
  #
  # Returns n x nc matrix of log-probabilities, or NULL on failure.
  # Training cells are hard-pinned: col = true_class gets 0, all others get -10000.
  compute_U_hat <- function(x_fit, y_fit, X_all) {
    y_factor <- factor(y_fit, levels = 1:nc)   # explicit levels -> columns in 1:nc order
    class_counts <- table(y_factor)

    # glmnet needs >= 1 sample per class; augment by duplication if needed
    if (any(class_counts < 1)) {
      warning("CALLR: some classes absent from LR training data")
      return(NULL)
    }

    fit <- tryCatch(
      cv.glmnet(x_fit, y_factor, family = "multinomial", alpha = 1),
      error = function(e) {
        warning("CALLR cv.glmnet failed: ", e$message)
        NULL
      }
    )
    if (is.null(fit)) return(NULL)

    pred <- predict(fit, X_all, s = "lambda.min", type = "response")
    if (is.array(pred)) pred <- pred[,,1] else pred <- as.matrix(pred)

    if (ncol(pred) != nc) {
      warning("CALLR: glmnet returned wrong number of class columns")
      return(NULL)
    }

    # Log-probabilities (U_hat in the original paper)
    pred    <- pmax(pred, 1e-10)
    U_hat   <- log(pred)

    # Hard-pin training cells: log-one-hot (-10000 for wrong class, 0 for true class)
    U_hat[T[,1], ] <- -10000
    for (k in 1:nrow(T)) {
      U_hat[T[k,1], T[k,2]] <- 0
    }

    return(U_hat)
  }

  # --- Initialise U (one-hot for training cells, uniform for unlabelled cells) ---
  U <- matrix(1 / nc, n, nc)
  U <- pin_training(U)

  # --- Initial logistic regression on training cells only ---
  U_hat <- compute_U_hat(X[T[,1], , drop = FALSE], T[,2], X)

  if (is.null(U_hat)) {
    # Fallback: pure label propagation with no LR guidance
    U_hat <- matrix(-10000 / nc, n, nc)
    U_hat[T[,1], ] <- -10000
    for (k in 1:nrow(T)) U_hat[T[k,1], T[k,2]] <- 0
  }

  # --- MBO parameters (from original CALLR paper) ---
  dt        <- 0.1   # step size
  NS        <- 3     # sub-steps per inner iteration
  max_outer <- 5     # outer iterations (refit LR on all cells)
  max_inner <- 100   # inner MBO iterations

  # ============================================================
  # Outer loop: binarise -> refit LR on all cells -> re-pin
  # ============================================================
  for (outer_iter in 1:max_outer) {
    U_outer_old <- U

    # ----------------------------------------------------------
    # Inner loop: MBO gradient steps with simplex projection
    # Original update: U <- U + (dt/NS) * (-L %*% U + u * U_hat)
    # ----------------------------------------------------------
    for (inner_iter in 1:max_inner) {
      U_inner_old <- U

      for (s in 1:NS) {
        U <- U + (dt / NS) * (-L %*% U + u * U_hat)
        # Project each row to probability simplex
        U <- t(apply(U, 1, ps))
      }

      # Re-pin training cells after each inner iteration
      U <- pin_training(U)

      if (max(abs(U - U_inner_old)) < 1e-5) break
    }

    # ----------------------------------------------------------
    # Binarise U: argmax -> one-hot; force training cells to true labels
    # ----------------------------------------------------------
    best_class            <- apply(U, 1, which.max)
    best_class[T[,1]]     <- T[,2]    # override with true labels

    U_bin <- matrix(0, n, nc)
    for (i in 1:n) U_bin[i, best_class[i]] <- 1

    # ----------------------------------------------------------
    # Refit LR on ALL cells using binarised labels
    # ----------------------------------------------------------
    U_hat_new <- compute_U_hat(X, best_class, X)
    if (!is.null(U_hat_new)) {
      U_hat <- U_hat_new
    }

    # Set U to binarised labels for next outer iteration
    U <- U_bin

    if (max(abs(U - U_outer_old)) < 1e-4) break
  }

  # --- Final predictions (argmax of U, mapped back to original label values) ---
  pred_indices <- apply(U, 1, which.max)   # indices into 1:nc
  return(pred_indices)                     # caller maps via label_vals[pred_indices]
}

#' Data preprocessing function
#' Preprocesses expression matrix (genes x cells) by removing zero-sum genes
#' and normalizing each cell column by its geometric mean of non-zero values.
preprocess <- function(X){
  # Remove genes with zero expression across all cells
  zs = which(apply(X,1,sum)==0)
  if(length(zs) > 0) {
    X = X[-zs,]
  }

  m = dim(X)[2]
  for(i in 1:m){
    a = X[,i]
    non_zero = a[a > 0]
    if(length(non_zero) > 0) {
      X[,i] = X[,i] / exp(mean(log(non_zero)))
    }
  }
  X = as.matrix(X)
  return(X)
}

#' Representative cell selection function
#'
#' Selects representative cells based on marker genes using TF-IDF scoring.
#' Enforces an exclusivity constraint: a cell is assigned to the cell type
#' for which it has the HIGHEST marker score (ties -> excluded).
#'
#' @param marker_file  Matrix, rows = cell types, cols = marker gene names (empty strings padded)
#' @param X            Matrix, genes x cells (output of preprocess())
#' @param cutoff       Numeric (0-1), global quantile threshold for TF-IDF zeroing (default 0.6)
#' @return List with training_set (2-col matrix: cell index, class label 1:nc) and label (cell type names)
representative <- function(marker_file, X, cutoff){
  n = dim(X)[1]   # genes
  m = dim(X)[2]   # cells

  # TF-IDF-like transformation
  rs = as.matrix(rowSums(X))   # gene totals
  cs = as.matrix(colSums(X))   # cell totals
  rs[rs == 0] = 1
  cs[cs == 0] = 1

  Xtf = X * as.vector(log(1 + m / rs)) / matrix(as.vector(cs), n, m, byrow = TRUE)

  # Zero out global bottom-quantile values
  cutoff_threshold = quantile(Xtf, cutoff, na.rm = TRUE)
  Xtf[Xtf < cutoff_threshold] = 0

  cell_types = rownames(marker_file)
  nc_rep     = length(cell_types)

  # --- Compute per-cell-type marker scores for ALL cells simultaneously ---
  all_scores = matrix(0, nc_rep, m)   # nc_rep x n_cells

  for (i in 1:nc_rep) {
    markers = marker_file[i, ]
    markers = markers[markers != "" & !is.na(markers)]
    if (length(markers) > 0) {
      marker_indices = which(rownames(X) %in% markers)
      if (length(marker_indices) > 0) {
        all_scores[i, ] = colSums(Xtf[marker_indices, , drop = FALSE])
      }
    }
  }

  # --- Exclusivity constraint: assign each cell to its argmax type ---
  # A cell is only a candidate representative for the type where it scores highest.
  best_type_per_cell = apply(all_scores, 2, which.max)   # length m

  representative_cells = c()
  labels = c()

  for (i in 1:nc_rep) {
    # Candidate cells: max score is for this type AND score is non-trivial
    candidates = which(best_type_per_cell == i & all_scores[i, ] > 0)

    if (length(candidates) < 4) {
      # Fallback 1: relax exclusivity — use top-scoring cells for this type
      relaxed = order(all_scores[i, ], decreasing = TRUE)[1:min(8, m)]
      candidates = relaxed[all_scores[i, relaxed] > 0]
    }

    if (length(candidates) < 4) {
      # Fallback 2: take top 8 cells regardless of score
      candidates = order(all_scores[i, ], decreasing = TRUE)[1:min(8, m)]
    }

    n_select = min(12, max(8, length(candidates)))
    selected = sample(candidates, n_select, replace = (length(candidates) < n_select))

    representative_cells = c(representative_cells, selected)
    labels               = c(labels, rep(i, n_select))
  }

  # Build training set matrix (cell indices, class labels)
  if (length(representative_cells) > 0) {
    training_set = cbind(representative_cells, labels)
  } else {
    warning("CALLR representative(): no representative cells found, using random selection")
    n_per_type   = 5
    training_set = matrix(0, nc_rep * n_per_type, 2)
    for (i in 1:nc_rep) {
      idx_start = (i - 1) * n_per_type + 1
      idx_end   = i * n_per_type
      training_set[idx_start:idx_end, 1] = sample(seq_len(m), n_per_type, replace = TRUE)
      training_set[idx_start:idx_end, 2] = i
    }
  }

  return(list(
    "training_set" = training_set,
    "label"        = cell_types
  ))
}

cat("CALLR core functions loaded successfully!\n")
