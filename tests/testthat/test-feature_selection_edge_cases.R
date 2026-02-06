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
