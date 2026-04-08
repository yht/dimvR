# Edge cases for feature selection

library(testthat)

test_that("single feature returns empty or minimal selection", {
  X <- data.frame(x1 = rnorm(20))
  res <- select_features_adaptive(X, target_var = "x1", method = "adaptive", min_features = 1)
  expect_true(is.list(res))
  expect_true(length(res$selected_features) <= 1)
})

test_that("all-missing target column handled gracefully", {
  X <- data.frame(x1 = rep(NA_real_, 30), x2 = rnorm(30), x3 = rnorm(30))
  res <- select_features_adaptive(X, target_var = "x1", method = "adaptive", min_features = 1)
  expect_true(is.list(res))
  expect_true(length(res$selected_features) >= 1)
})

test_that("very small n does not error", {
  X <- data.frame(x1 = rnorm(2), x2 = rnorm(2), x3 = rnorm(2))
  res <- select_features_adaptive(X, target_var = "x1", method = "mi", min_features = 1)
  expect_true(is.list(res))
})

test_that("invalid feature-selection arguments fail fast", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))

  expect_error(
    select_features_adaptive(X, target_var = "x1", min_features = 3, max_features = 2),
    "min_features must be <= max_features"
  )

  expect_error(
    select_features_adaptive(X, target_var = "x1", method = "fixed", threshold = 2),
    "threshold must be between 0 and 1"
  )

  expect_error(
    select_features_adaptive(X, target_var = "x1", method = "mi", nbins = 1),
    "nbins must be >= 2"
  )
})

test_that("min_features larger than candidate count is handled safely", {
  X <- data.frame(x1 = rnorm(25), x2 = rnorm(25), x3 = rnorm(25))
  res <- select_features_adaptive(
    X,
    target_var = "x1",
    method = "adaptive",
    min_features = 10,
    max_features = 10
  )

  expect_true(is.list(res))
  expect_equal(res$candidate_count, 2L)
  expect_lte(length(res$selected_features), 2)
})
