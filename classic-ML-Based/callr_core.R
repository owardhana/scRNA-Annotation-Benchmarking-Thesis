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
#' Projects a vector to a simplex
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
#' Performs semi-supervised classification using Laplacian and logistic regression
callr <- function(X, u, T){

  n = dim(X)[1]
  m = dim(X)[2]
  nc = length(unique(T[,2]))

  # Initialize F matrix
  F = matrix(0, n, nc)
  for(i in 1:nc){
    F[T[which(T[,2]==sort(unique(T[,2]))[i]),1],i] = 1
  }

  # Compute distance and weight matrix with improved sigma calculation
  D = dist2(X)
  # Use median of k nearest neighbors for more robust sigma calculation
  k_neighbors = min(10, floor(n/4))  # Use up to 10 neighbors or n/4, whichever is smaller
  sigma = apply(D, 1, function(row) {
    # Get k-th nearest neighbor distance, avoiding exact zeros
    sorted_dists = sort(row[row > 0])
    if(length(sorted_dists) >= k_neighbors) {
      return(sorted_dists[k_neighbors])
    } else if(length(sorted_dists) > 0) {
      return(sorted_dists[length(sorted_dists)])
    } else {
      return(1)  # fallback value
    }
  })

  # Ensure sigma values are reasonable (not too small)
  sigma = pmax(sigma, quantile(sigma, 0.1, na.rm = TRUE))

  W = exp(-D/(sigma%*%t(sigma)))
  diag(W) = 0

  # Compute Laplacian
  d = apply(W,1,sum)
  L = diag(d) - W

  # Initial logistic regression
  x.train = X[T[,1],]
  y.train = T[,2]

  # Labels will be converted to factors in the validation section above

  # Validate and augment training set to meet glmnet requirements
  class_counts <- table(y.train)
  min_samples_needed <- 8

  # Check if any class has insufficient samples
  insufficient_classes <- names(class_counts)[class_counts < min_samples_needed]

  if(length(insufficient_classes) > 0) {
    # Augment training set by duplicating samples from insufficient classes
    augmented_indices <- T[,1]
    augmented_labels <- T[,2]

    for(cls_idx in insufficient_classes) {
      cls_samples <- which(T[,2] == as.numeric(cls_idx))
      current_count <- length(cls_samples)
      needed_samples <- min_samples_needed - current_count

      if(needed_samples > 0 && current_count > 0) {
        # Duplicate existing samples with small noise to reach minimum
        duplicates_needed <- ceiling(needed_samples / current_count)
        for(dup in 1:duplicates_needed) {
          if(length(augmented_indices) + current_count <= n * 2) {  # Safety limit
            augmented_indices <- c(augmented_indices, T[cls_samples, 1])
            augmented_labels <- c(augmented_labels, T[cls_samples, 2])
          }
          if(length(augmented_labels[augmented_labels == as.numeric(cls_idx)]) >= min_samples_needed) break
        }
      }
    }

    # Update training data with augmented set
    x.train = X[augmented_indices,]
    y.train = as.factor(augmented_labels)
  } else {
    # Original training data is sufficient
    x.train = X[T[,1],]
    y.train = as.factor(T[,2])
  }

  # Fit multinomial logistic regression with balanced classes
  # Create class weights to handle imbalanced training data
  class_counts <- table(y.train)

  # Only proceed if we have adequate samples per class
  if(min(class_counts) >= 3) {  # Minimum viable threshold
    class_weights <- rep(1, length(class_counts))
    names(class_weights) <- names(class_counts)

    # Calculate inverse frequency weights for balancing
    total_samples <- sum(class_counts)
    for(cls in names(class_counts)) {
      class_weights[cls] <- total_samples / (length(class_counts) * class_counts[cls])
    }

    # Apply weights to training samples
    sample_weights <- class_weights[as.character(y.train)]

    fit = cv.glmnet(x.train, y.train, family="multinomial", alpha=1, weights=sample_weights)
  } else {
    # Fallback: use simple glmnet without weights if still insufficient
    warning("Insufficient training samples for some classes, using unweighted glmnet")
    fit = cv.glmnet(x.train, y.train, family="multinomial", alpha=1)
  }

  # Predict on all data
  pred = predict(fit, X, s="lambda.min", type="response")

  # Initialize F with predictions
  if(is.array(pred)) {
    F = pred[,,1]
  } else {
    F = as.matrix(pred)
  }

  # Ensure F has correct dimensions
  if(ncol(F) != nc) {
    # If predictions don't match expected classes, initialize with uniform
    F = matrix(1/nc, n, nc)
    for(i in 1:nc){
      F[T[which(T[,2]==sort(unique(T[,2]))[i]),1],i] = 1
    }
  }

  # Iterative optimization
  for(iter in 1:100) {
    F_old = F

    # Update F using graph regularization
    for(i in 1:n) {
      if(!(i %in% T[,1])) {  # Only update unlabeled points
        # Graph-based update
        numerator = W[i,] %*% F
        denominator = sum(W[i,])
        if(denominator > 0) {
          F[i,] = numerator / denominator
        }
        # Project to simplex
        F[i,] = ps(F[i,] + u * (L[i,] %*% F))
      }
    }

    # Check convergence
    if(max(abs(F - F_old)) < 1e-6) {
      break
    }
  }

  # Return predicted labels with balanced selection
  # Instead of just taking max, use a more balanced approach
  predictions <- apply(F, 1, function(row) {
    # If prediction is very confident (>0.7), use it
    max_prob <- max(row)
    if(max_prob > 0.7) {
      return(which.max(row))
    } else {
      # For uncertain predictions, consider balance
      # Add small random noise to break ties more fairly
      noisy_row <- row + runif(length(row), -0.01, 0.01)
      return(which.max(noisy_row))
    }
  })

  return(predictions)
}

#' Data preprocessing function
#' Preprocesses expression matrix by removing zero genes and normalizing
preprocess <- function(X){
  # Remove genes with zero expression across all cells
  zs = which(apply(X,1,sum)==0)
  if(length(zs) > 0) {
    X = X[-zs,]
  }

  m = dim(X)[2]
  for(i in 1:m){
    a = X[,i]
    a = as.matrix(a)
    non_zero = a[a > 0]
    if(length(non_zero) > 0) {
      X[,i] = X[,i]/exp(mean(log(non_zero)))
    }
  }
  X = as.matrix(X)
  return(X)
}

#' Representative cell selection function
#' Selects representative cells based on marker genes
representative <- function(marker_file, X, cutoff){
  n = dim(X)[1]
  m = dim(X)[2]

  # Calculate TF-IDF transformation matrix
  rs = as.matrix(apply(X,1,sum))
  cs = as.matrix(apply(X,2,sum))

  # Avoid division by zero
  rs[rs == 0] = 1
  cs[cs == 0] = 1

  T = X * (log(1 + m/rs) %*% t(1/cs))

  # Calculate cutoff threshold
  cutoff_threshold = quantile(T, cutoff, na.rm = TRUE)
  T[T < cutoff_threshold] = 0

  # Find representative cells for each cell type
  cell_types = rownames(marker_file)
  representative_cells = c()
  labels = c()

  for(i in 1:length(cell_types)) {
    markers = marker_file[i,]
    markers = markers[markers != "" & !is.na(markers)]

    if(length(markers) > 0) {
      # Find genes that match markers
      marker_indices = which(rownames(X) %in% markers)

      if(length(marker_indices) > 0) {
        # Calculate score for each cell based on marker expression
        marker_scores = apply(T[marker_indices, , drop=FALSE], 2, sum)

        # Select top cells above cutoff
        top_cells = which(marker_scores >= quantile(marker_scores, cutoff, na.rm = TRUE))

        if(length(top_cells) > 0) {
          # Select 8-12 representative cells per type to meet glmnet requirements
          n_select = min(12, max(8, length(top_cells)))

          if(length(top_cells) < 8) {
            # If very few cells meet cutoff, progressively lower the threshold
            relaxed_cutoff = max(0.2, cutoff - 0.4)
            top_cells = which(marker_scores >= quantile(marker_scores, relaxed_cutoff, na.rm = TRUE))

            # If still insufficient, take top scoring cells regardless of cutoff
            if(length(top_cells) < 8) {
              top_cells = order(marker_scores, decreasing = TRUE)[1:min(8, m)]
            }
            n_select = min(10, max(8, length(top_cells)))
          }
          selected = sample(top_cells, n_select)

          representative_cells = c(representative_cells, selected)
          labels = c(labels, rep(i, n_select))
        }
      }
    }
  }

  # Create training set matrix
  if(length(representative_cells) > 0) {
    training_set = cbind(representative_cells, labels)
  } else {
    # If no representatives found, select random cells
    warning("No representative cells found, selecting random cells")
    n_per_type = 5
    training_set = matrix(0, length(cell_types) * n_per_type, 2)
    for(i in 1:length(cell_types)) {
      idx_start = (i-1) * n_per_type + 1
      idx_end = i * n_per_type
      training_set[idx_start:idx_end, 1] = sample(1:m, n_per_type)
      training_set[idx_start:idx_end, 2] = i
    }
  }

  return(list(
    "training_set" = training_set,
    "label" = cell_types
  ))
}

cat("CALLR core functions loaded successfully!\n")