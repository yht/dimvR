# Test suite for feature selection module
# Phase 1, Days 1-3

library(testthat)

test_that("select_features_adaptive works with simple data", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50),
    x4 = rnorm(50)
  )
  
  # Make x2 highly correlated with x1
  X$x2 <- X$x1 + rnorm(50, sd = 0.1)
  
  result <- select_features_adaptive(X, target_var = 1, method = "adaptive")
  
  expect_true(is.list(result))
  expect_true("selected_features" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true(length(result$selected_features) >= 2)
  expect_true(2 %in% result$selected_features)  # x2 should be selected
})


test_that("select_features_adaptive respects min_features", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  
  result <- select_features_adaptive(X, target_var = 1, 
                                     method = "adaptive", min_features = 2)
  
  expect_true(length(result$selected_features) >= 2)
})


test_that("select_features_adaptive respects max_features", {
  set.seed(123)
  X <- data.frame(matrix(rnorm(50 * 10), ncol = 10))
  
  result <- select_features_adaptive(X, target_var = 1, 
                                     method = "adaptive", max_features = 3)
  
  expect_true(length(result$selected_features) <= 3)
})


test_that("fixed threshold selection works", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  X$x2 <- 0.8 * X$x1 + rnorm(50, sd = 0.3)
  
  result <- select_features_adaptive(X, target_var = 1, 
                                     method = "fixed", threshold = 0.5)
  
  expect_true(2 %in% result$selected_features)
})


test_that("mutual information selection works", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    x3 = rnorm(100)
  )
  
  result <- select_features_adaptive(X, target_var = 1, method = "mi")
  
  expect_true(is.list(result))
  expect_true(length(result$selected_features) >= 2)
})


test_that("hybrid selection combines correlation and MI", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    x3 = rnorm(100)
  )
  
  result <- select_features_adaptive(X, target_var = 1, method = "hybrid")
  
  expect_true(is.list(result))
  expect_true("method_info" %in% names(result))
  expect_true("correlation_scores" %in% names(result$method_info))
  expect_true("mi_scores" %in% names(result$method_info))
})


test_that("feature selection handles missing data", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  
  # Introduce missing values
  X$x2[1:5] <- NA
  X$x3[6:10] <- NA
  
  result <- select_features_adaptive(X, target_var = 1, method = "adaptive")
  
  expect_true(is.list(result))
  expect_true(length(result$selected_features) >= 1)
})


test_that("feature selection works with character target_var", {
  set.seed(123)
  X <- data.frame(
    a = rnorm(50),
    b = rnorm(50),
    c = rnorm(50)
  )
  
  result <- select_features_adaptive(X, target_var = "a", method = "adaptive")
  
  expect_true(is.list(result))
  expect_equal(result$target_var, 1)
})


test_that("compute_simple_mi returns non-negative values", {
  set.seed(123)
  x <- rnorm(100)
  y <- rnorm(100)
  
  mi <- compute_simple_mi(x, y, nbins = 5)
  
  expect_true(mi >= 0)
  expect_true(is.finite(mi))
})


test_that("compute_simple_mi is higher for correlated variables", {
  set.seed(123)
  x <- rnorm(100)
  y_random <- rnorm(100)
  y_correlated <- x + rnorm(100, sd = 0.1)
  
  mi_random <- compute_simple_mi(x, y_random)
  mi_correlated <- compute_simple_mi(x, y_correlated)
  
  expect_true(mi_correlated > mi_random)
})

test_that("select_features_adaptive returns stable metadata contract", {
  set.seed(123)
  X <- data.frame(
    a = rnorm(60),
    b = rnorm(60),
    c = rnorm(60),
    d = rnorm(60)
  )

  res <- select_features_adaptive(
    X,
    target_var = "a",
    method = "hybrid",
    min_features = 2,
    max_features = 3
  )

  expect_true("target_index" %in% names(res))
  expect_true("target_name" %in% names(res))
  expect_true("method" %in% names(res))
  expect_true("n_selected" %in% names(res))
  expect_true("candidate_count" %in% names(res))
  expect_true("fallback_used" %in% names(res))
  expect_true("ranking" %in% names(res))

  expect_equal(res$target_index, 1L)
  expect_equal(res$target_name, "a")
  expect_equal(res$method, "hybrid")
  expect_equal(res$candidate_count, 3L)
  expect_equal(res$n_selected, length(res$selected_features))
  expect_equal(length(res$selected_names), length(res$selected_features))
  expect_true(is.logical(res$fallback_used) && length(res$fallback_used) == 1)
  expect_s3_class(res$ranking, "data.frame")
  expect_true(all(c("feature", "final_score", "rank") %in% names(res$ranking)))
  expect_equal(nrow(res$ranking), res$candidate_count)
})
