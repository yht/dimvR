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

#' Adaptive Feature Selection for DIMV (Correlation-based)
#'
#' Selects a subset of predictor features for imputing a target variable,
#' using a simple adaptive correlation-based threshold. The goal is to
#' reduce noise and multicollinearity in the conditional models while
#' maintaining minimal intrusion into the DIMV framework.
#'
#' @keywords internal
#' @importFrom stats cor median
select_features_adaptive <- function(X,
                                     target_var,
                                     min_features = 2,
                                     max_features = 5,
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
  
  candidates <- setdiff(seq_len(p), target_var)
  target_col <- X[, target_var]
  
  # Pairwise correlations with target
  cors <- sapply(candidates, function(j) {
    complete_both <- !is.na(target_col) & !is.na(X[, j])
    if (sum(complete_both) < 3) return(0)  # Not enough data for correlation
    cor(target_col[complete_both], X[complete_both, j], use = "complete.obs")
  })
  
  abs_cors <- abs(cors)
  names(abs_cors) <- colnames(X)[candidates]
  
  # Adaptive threshold: median absolute correlation among non-zero
  non_zero <- abs_cors[abs_cors > 0]
  median_cor <- if (length(non_zero)) stats::median(non_zero) else 0
  
  sorted <- sort(abs_cors, decreasing = TRUE)
  
  # Rule: use median but ensure >= min_features
  if (sum(abs_cors >= median_cor) < min_features) {
    threshold_used <- sorted[min_features]
  } else {
    threshold_used <- median_cor
  }
  
  # Cap by max_features
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_features) {
    selected <- names(sorted[seq_len(max_features)])
  }
  
  if (verbose) {
    cat("Adaptive FS: threshold =", round(threshold_used, 3),
        "| selected", length(selected), "features\n")
  }
  
  list(
    selected_features = selected,
    scores = abs_cors,
    threshold_used = threshold_used,
    method_info = list(
      strategy = "median_adaptive",
      median_correlation = median_cor
    )
  )
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

