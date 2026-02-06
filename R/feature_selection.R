# ===========================================================
# File: R/feature_selection.R
# Intelligent Feature Selection for DIMV
# Status: EXPERIMENTAL (Phase 1 - Internal API)
#
# This module is currently experimental and not part of the
# official stable DIMV API. Interfaces may change in future
# releases. DO NOT rely on these functions for production use.
#
# To promote functions to Public API in the future:
#   - Add @export tags
#   - Add unit tests
#   - Add vignette examples
# ===========================================================

# Internal helper to show experimental warning (once per session)
.dimv_fs_warn_once <- local({
  warned <- FALSE
  function() {
    if (!warned) {
      warning("[dimvR] Feature Selection module is EXPERIMENTAL and may change without notice.",
              call. = FALSE)
      warned <<- TRUE
    }
  }
})

#' Adaptive Feature Selection for DIMV (Correlation/MI/Hybrid)
#'
#' Selects a subset of predictor features for imputing a target variable.
#' Supports adaptive (correlation), fixed-threshold (correlation),
#' mutual information, and hybrid (correlation + MI) methods.
#'
#' @keywords internal
#' @importFrom stats cor median
select_features_adaptive <- function(X,
                                     target_var,
                                     min_features = 2,
                                     max_features = 5,
                                     method = c("adaptive", "fixed", "mi", "hybrid"),
                                     threshold = NULL,
                                     nbins = 5,
                                     verbose = FALSE) {
  
  .dimv_fs_warn_once()
  
  X <- as.data.frame(lapply(X, as.numeric))
  p <- ncol(X)
  
  # Convert target_var to index if name provided
  if (is.character(target_var)) {
    target_var <- which(colnames(X) == target_var)
  }
  
  if (length(target_var) != 1 || target_var < 1 || target_var > p) {
    stop("target_var must be a valid column index or name.")
  }
  
  method <- match.arg(method)
  candidates <- setdiff(seq_len(p), target_var)
  target_col <- X[, target_var]
  
  if (length(candidates) == 0) {
    return(list(
      selected_features = integer(0),
      selected_names = character(0),
      scores = numeric(0),
      target_var = target_var,
      threshold_used = NA_real_,
      method_info = list(
        method = match.arg(method),
        correlation_scores = numeric(0),
        mi_scores = NULL
      )
    ))
  }
  
  # Pairwise correlations with target
  cors <- vapply(candidates, function(j) {
    complete_both <- !is.na(target_col) & !is.na(X[, j])
    if (sum(complete_both) < 3) return(0)  # Not enough data for correlation
    cor(target_col[complete_both], X[complete_both, j], use = "complete.obs")
  }, numeric(1))
  
  abs_cors <- abs(cors)
  names(abs_cors) <- colnames(X)[candidates]
  
  mi_scores <- NULL
  selected_names <- character(0)
  threshold_used <- NA_real_
  
  if (method == "adaptive") {
    res <- select_adaptive_threshold(abs_cors,
                                     min_features = min_features,
                                     max_features = max_features)
    selected_names <- res$selected
    threshold_used <- res$threshold_used
  } else if (method == "fixed") {
    if (is.null(threshold)) stop("threshold must be provided for method = 'fixed'.")
    res <- select_fixed_threshold(abs_cors,
                                  threshold = threshold,
                                  min_features = min_features,
                                  max_features = max_features)
    selected_names <- res$selected
    threshold_used <- res$threshold_used
  } else if (method == "mi") {
    res <- select_mutual_information(X,
                                     target_var = target_var,
                                     min_features = min_features,
                                     max_features = max_features,
                                     nbins = nbins)
    selected_names <- res$selected
    mi_scores <- res$mi_scores
    threshold_used <- res$threshold_used
  } else if (method == "hybrid") {
    res <- select_hybrid(X,
                         target_var = target_var,
                         min_features = min_features,
                         max_features = max_features,
                         nbins = nbins)
    selected_names <- res$selected
    mi_scores <- res$mi_scores
    threshold_used <- res$threshold_used
  }
  
  if (length(selected_names) == 0) {
    selected_names <- names(sort(abs_cors, decreasing = TRUE))[seq_len(min_features)]
  }
  
  if (verbose) {
    cat("Feature selection:", method,
        "| selected", length(selected_names), "features\n")
  }
  
  selected_idx <- match(selected_names, colnames(X))
  selected_idx <- selected_idx[!is.na(selected_idx)]
  
  list(
    selected_features = selected_idx,
    selected_names = selected_names,
    scores = abs_cors,
    target_var = target_var,
    threshold_used = threshold_used,
    method_info = list(
      method = method,
      correlation_scores = abs_cors,
      mi_scores = mi_scores
    )
  )
}

#' Adaptive threshold selection (correlation-based)
#' @keywords internal
select_adaptive_threshold <- function(abs_cors,
                                      min_features = 2,
                                      max_features = 5) {
  non_zero <- abs_cors[abs_cors > 0]
  median_cor <- if (length(non_zero)) stats::median(non_zero) else 0
  sorted <- sort(abs_cors, decreasing = TRUE)
  
  if (length(sorted) == 0) {
    return(list(selected = character(0), threshold_used = NA_real_))
  }
  
  if (sum(abs_cors >= median_cor) < min_features) {
    threshold_used <- sorted[min_features]
  } else {
    threshold_used <- median_cor
  }
  
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_features) {
    selected <- names(sorted[seq_len(max_features)])
  }
  
  list(selected = selected, threshold_used = threshold_used)
}

#' Fixed threshold selection (correlation-based)
#' @keywords internal
select_fixed_threshold <- function(abs_cors,
                                   threshold = 0.3,
                                   min_features = 2,
                                   max_features = 5) {
  sorted <- sort(abs_cors, decreasing = TRUE)
  
  if (length(sorted) == 0) {
    return(list(selected = character(0), threshold_used = threshold))
  }
  
  selected <- names(sorted[sorted >= threshold])
  if (length(selected) < min_features) {
    selected <- names(sorted[seq_len(min_features)])
  }
  if (length(selected) > max_features) {
    selected <- names(sorted[seq_len(max_features)])
  }
  
  list(selected = selected, threshold_used = threshold)
}

#' Mutual information selection
#' @keywords internal
select_mutual_information <- function(X,
                                      target_var,
                                      min_features = 2,
                                      max_features = 5,
                                      nbins = 5) {
  X <- as.data.frame(lapply(X, as.numeric))
  p <- ncol(X)
  candidates <- setdiff(seq_len(p), target_var)
  target_col <- X[, target_var]
  
  mi <- sapply(candidates, function(j) {
    complete_both <- !is.na(target_col) & !is.na(X[, j])
    if (sum(complete_both) < 3) return(0)
    compute_simple_mi(target_col[complete_both], X[complete_both, j], nbins = nbins)
  })
  
  mi_scores <- mi
  names(mi_scores) <- colnames(X)[candidates]
  sorted <- sort(mi_scores, decreasing = TRUE)
  
  if (length(sorted) == 0) {
    return(list(selected = character(0), threshold_used = NA_real_, mi_scores = mi_scores))
  }
  
  threshold_used <- sorted[min_features]
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_features) {
    selected <- names(sorted[seq_len(max_features)])
  }
  
  list(selected = selected, threshold_used = threshold_used, mi_scores = mi_scores)
}

#' Hybrid selection (correlation + MI)
#' @keywords internal
select_hybrid <- function(X,
                          target_var,
                          min_features = 2,
                          max_features = 5,
                          nbins = 5) {
  X <- as.data.frame(lapply(X, as.numeric))
  p <- ncol(X)
  candidates <- setdiff(seq_len(p), target_var)
  target_col <- X[, target_var]
  
  cors <- sapply(candidates, function(j) {
    complete_both <- !is.na(target_col) & !is.na(X[, j])
    if (sum(complete_both) < 3) return(0)
    stats::cor(target_col[complete_both], X[complete_both, j], use = "complete.obs")
  })
  abs_cors <- abs(cors)
  names(abs_cors) <- colnames(X)[candidates]
  
  mi_scores <- sapply(candidates, function(j) {
    complete_both <- !is.na(target_col) & !is.na(X[, j])
    if (sum(complete_both) < 3) return(0)
    compute_simple_mi(target_col[complete_both], X[complete_both, j], nbins = nbins)
  })
  names(mi_scores) <- colnames(X)[candidates]
  
  # Normalize to [0,1]
  norm <- function(v) {
    rng <- range(v, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] - rng[1] == 0) {
      return(rep(0, length(v)))
    }
    (v - rng[1]) / (rng[2] - rng[1])
  }
  
  score <- 0.5 * norm(abs_cors) + 0.5 * norm(mi_scores)
  names(score) <- colnames(X)[candidates]
  sorted <- sort(score, decreasing = TRUE)
  
  if (length(sorted) == 0) {
    return(list(selected = character(0), threshold_used = NA_real_, mi_scores = mi_scores))
  }
  
  threshold_used <- sorted[min_features]
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_features) {
    selected <- names(sorted[seq_len(max_features)])
  }
  
  list(selected = selected, threshold_used = threshold_used, mi_scores = mi_scores)
}


#' Simple Mutual Information Estimator (Experimental)
#'
#' Computes MI using histogram-based binning. This is a basic implementation
#' suitable for continuous variables. Not yet optimized and may be replaced by
#' a kernel-based method in future releases.
#'
#' @keywords internal
compute_simple_mi <- function(x, y, nbins = 5) {
  
  .dimv_fs_warn_once()
  
  # Discretize into bins
  x_bins <- cut(x, breaks = nbins, labels = FALSE, include.lowest = TRUE)
  y_bins <- cut(y, breaks = nbins, labels = FALSE, include.lowest = TRUE)
  
  # Joint distribution
  joint_table <- table(x_bins, y_bins)
  p_xy <- joint_table / sum(joint_table)
  
  # Marginals
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  
  # MI = sum p(x,y) * log(p(x,y) / (p(x) * p(y)))
  mi <- 0
  for (i in seq_along(p_x)) {
    for (j in seq_along(p_y)) {
      if (p_xy[i, j] > 0) {
        mi <- mi + p_xy[i, j] * log(p_xy[i, j] / (p_x[i] * p_y[j]))
      }
    }
  }
  
  return(max(0, mi))  # MI should be non-negative
}

