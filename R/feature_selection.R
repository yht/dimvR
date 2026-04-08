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
#' @param X A numeric data frame or matrix.
#' @param target_var Target variable (name or index).
#' @param min_features Minimum number of features to select.
#' @param max_features Maximum number of features to select.
#' @param method Selection method: adaptive, fixed, mi, or hybrid.
#' @param threshold Fixed correlation threshold (for method = "fixed").
#' @param nbins Number of bins for MI estimation (methods "mi" and "hybrid").
#' @param verbose If TRUE, prints selection summary.
#'
#' @export
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

  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data.frame or matrix.")
  }

  X <- as.data.frame(lapply(X, as.numeric))
  if (ncol(X) == 0) {
    stop("X must have at least one column.")
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }

  p <- ncol(X)

  if (!is.numeric(min_features) || length(min_features) != 1 || is.na(min_features)) {
    stop("min_features must be a single numeric value.")
  }
  if (!is.numeric(max_features) || length(max_features) != 1 || is.na(max_features)) {
    stop("max_features must be a single numeric value.")
  }
  min_features <- as.integer(min_features)
  max_features <- as.integer(max_features)

  if (min_features < 1) {
    stop("min_features must be >= 1.")
  }
  if (max_features < 1) {
    stop("max_features must be >= 1.")
  }
  if (min_features > max_features) {
    stop("min_features must be <= max_features.")
  }

  method <- match.arg(method)

  if (method == "fixed") {
    if (is.null(threshold) || !is.numeric(threshold) || length(threshold) != 1 || is.na(threshold)) {
      stop("threshold must be a single numeric value for method = 'fixed'.")
    }
    if (threshold < 0 || threshold > 1) {
      stop("threshold must be between 0 and 1 for method = 'fixed'.")
    }
  }

  if (method %in% c("mi", "hybrid")) {
    if (!is.numeric(nbins) || length(nbins) != 1 || is.na(nbins)) {
      stop("nbins must be a single numeric value for method 'mi' and 'hybrid'.")
    }
    nbins <- as.integer(nbins)
    if (nbins < 2) {
      stop("nbins must be >= 2 for method 'mi' and 'hybrid'.")
    }
  }

  # Convert target_var to index if name provided
  if (is.character(target_var)) {
    target_var <- which(colnames(X) == target_var)
  }

  if (length(target_var) != 1 || target_var < 1 || target_var > p) {
    stop("target_var must be a valid column index or name.")
  }
  target_var <- as.integer(target_var)

  candidates <- setdiff(seq_len(p), target_var)
  target_col <- X[, target_var]

  if (length(candidates) == 0) {
    return(list(
      selected_features = integer(0),
      selected_names = character(0),
      scores = numeric(0),
      target_var = target_var,
      target_index = target_var,
      target_name = colnames(X)[target_var],
      method = method,
      n_selected = 0L,
      candidate_count = 0L,
      threshold_used = NA_real_,
      fallback_used = FALSE,
      ranking = data.frame(
        feature = character(0),
        correlation_score = numeric(0),
        mi_score = numeric(0),
        hybrid_score = numeric(0),
        final_score = numeric(0),
        rank = integer(0),
        stringsAsFactors = FALSE
      ),
      method_info = list(
        method = method,
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
  fallback_used <- FALSE
  ranking_final <- abs_cors
  ranking_hybrid <- NULL

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
    ranking_final <- mi_scores
  } else if (method == "hybrid") {
    res <- select_hybrid(X,
                         target_var = target_var,
                         min_features = min_features,
                         max_features = max_features,
                         nbins = nbins)
    selected_names <- res$selected
    mi_scores <- res$mi_scores
    threshold_used <- res$threshold_used
    ranking_hybrid <- res$hybrid_scores
    ranking_final <- ranking_hybrid
  }

  if (length(selected_names) == 0) {
    fallback_n <- min(max(min_features, 1L), length(abs_cors))
    selected_names <- names(sort(abs_cors, decreasing = TRUE))[seq_len(fallback_n)]
    fallback_used <- TRUE
  }

  if (verbose) {
    cat("Feature selection:", method,
        "| selected", length(selected_names), "features\n")
  }

  selected_idx <- match(selected_names, colnames(X))
  selected_idx <- selected_idx[!is.na(selected_idx)]
  selected_names <- colnames(X)[selected_idx]
  selected_idx <- as.integer(selected_idx)

  candidate_names <- colnames(X)[candidates]
  correlation_for_ranking <- abs_cors[candidate_names]
  mi_for_ranking <- if (is.null(mi_scores)) rep(NA_real_, length(candidate_names)) else mi_scores[candidate_names]
  hybrid_for_ranking <- if (is.null(ranking_hybrid)) rep(NA_real_, length(candidate_names)) else ranking_hybrid[candidate_names]
  final_for_ranking <- ranking_final[candidate_names]

  ranking <- data.frame(
    feature = candidate_names,
    correlation_score = as.numeric(correlation_for_ranking),
    mi_score = as.numeric(mi_for_ranking),
    hybrid_score = as.numeric(hybrid_for_ranking),
    final_score = as.numeric(final_for_ranking),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(-ranking$final_score, ranking$feature), , drop = FALSE]
  ranking$rank <- seq_len(nrow(ranking))

  list(
    selected_features = selected_idx,
    selected_names = selected_names,
    scores = abs_cors,
    target_var = target_var,
    target_index = target_var,
    target_name = colnames(X)[target_var],
    method = method,
    n_selected = length(selected_idx),
    candidate_count = length(candidates),
    threshold_used = threshold_used,
    fallback_used = fallback_used,
    ranking = ranking,
    method_info = list(
      method = method,
      correlation_scores = abs_cors,
      mi_scores = mi_scores,
      hybrid_scores = ranking_hybrid,
      final_scores = ranking_final
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

  min_n <- min(max(min_features, 1L), length(sorted))
  max_n <- min(max(max_features, min_n), length(sorted))

  if (sum(abs_cors >= median_cor) < min_n) {
    threshold_used <- sorted[min_n]
  } else {
    threshold_used <- median_cor
  }

  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_n) {
    selected <- names(sorted[seq_len(max_n)])
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

  min_n <- min(max(min_features, 1L), length(sorted))
  max_n <- min(max(max_features, min_n), length(sorted))

  selected <- names(sorted[sorted >= threshold])
  if (length(selected) < min_n) {
    selected <- names(sorted[seq_len(min_n)])
  }
  if (length(selected) > max_n) {
    selected <- names(sorted[seq_len(max_n)])
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

  min_n <- min(max(min_features, 1L), length(sorted))
  max_n <- min(max(max_features, min_n), length(sorted))
  threshold_used <- sorted[min_n]
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_n) {
    selected <- names(sorted[seq_len(max_n)])
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
    return(list(
      selected = character(0),
      threshold_used = NA_real_,
      mi_scores = mi_scores,
      hybrid_scores = score
    ))
  }

  min_n <- min(max(min_features, 1L), length(sorted))
  max_n <- min(max(max_features, min_n), length(sorted))
  threshold_used <- sorted[min_n]
  selected <- names(sorted[sorted >= threshold_used])
  if (length(selected) > max_n) {
    selected <- names(sorted[seq_len(max_n)])
  }

  list(
    selected = selected,
    threshold_used = threshold_used,
    mi_scores = mi_scores,
    hybrid_scores = score
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

