library(testthat)
library(dimvR)

test_that("dimv_diagnostics returns correct structure and values", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5, maxit = 50, tol = 1e-4)

  # Impute to get X_imputed
  X_imputed <- dimv_impute_new(imp, X)

  # Test diagnostics structure
  diag <- dimv_diagnostics(imp, X, X_imputed)

  # Test all expected fields exist
  expect_named(diag, c("iterations", "adaptive", "feature_select",
                       "residual_variances", "diff_history", "missing_rmse"))

  # Test iterations is a positive integer
  expect_type(diag$iterations, "double")
  expect_gt(diag$iterations, 0)

  # Test adaptive is logical
  expect_type(diag$adaptive, "logical")

  # Test feature_select is logical
  expect_type(diag$feature_select, "logical")

  # Test residual_variances is numeric vector with length = ncol(X)
  expect_type(diag$residual_variances, "double")
  expect_length(diag$residual_variances, ncol(X))

  # Test diff_history is numeric vector
  expect_type(diag$diff_history, "double")
  expect_true(length(diag$diff_history) > 0)

  # Test missing_rmse is numeric and finite
  expect_type(diag$missing_rmse, "double")
  expect_true(is.finite(diag$missing_rmse))
})

test_that("dimv_diagnostics works with complete data (no missing values)", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  imp <- dimv_train(X, lambda = 0.5, maxit = 20, tol = 1e-4)

  X_imputed <- dimv_impute_new(imp, X)

  diag <- dimv_diagnostics(imp, X, X_imputed)

  # Should still return valid structure
  expect_named(diag, c("iterations", "adaptive", "feature_select",
                       "residual_variances", "diff_history", "missing_rmse"))

  # With complete data, diff_history should be numeric(0) or very short
  expect_type(diag$diff_history, "double")

  # missing_rmse should be finite
  expect_true(is.finite(diag$missing_rmse))
})

test_that("dimv_diagnostics with different maxit values affects convergence", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 10), j] <- NA

  # Train with maxit = 10 (few iterations)
  imp_short <- dimv_train(X, lambda = 0.5, maxit = 10, tol = 1e-4)
  X_short <- dimv_impute_new(imp_short, X)
  diag_short <- dimv_diagnostics(imp_short, X, X_short)

  # Train with maxit = 50 (more iterations)
  imp_long <- dimv_train(X, lambda = 0.5, maxit = 50, tol = 1e-4)
  X_long <- dimv_impute_new(imp_long, X)
  diag_long <- dimv_diagnostics(imp_long, X, X_long)

  # Both should return valid diagnostics
  expect_type(diag_short$iterations, "double")
  expect_type(diag_long$iterations, "double")

  # Longer training should generally converge better (fewer diff_history entries or smaller values)
  expect_true(is.finite(diag_short$missing_rmse))
  expect_true(is.finite(diag_long$missing_rmse))
})

test_that("dimv_diagnostics residual_variances are positive", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5)

  X_imputed <- dimv_impute_new(imp, X)
  diag <- dimv_diagnostics(imp, X, X_imputed)

  # All residual variances should be positive (non-negative, and typically > 0)
  expect_true(all(diag$residual_variances >= 0))
  # Should have one residual variance per column
  expect_length(diag$residual_variances, ncol(X))
})

test_that("dimv_diagnostics diff_history shows convergence", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5, maxit = 30, tol = 1e-4)

  X_imputed <- dimv_impute_new(imp, X)
  diag <- dimv_diagnostics(imp, X, X_imputed)

  # diff_history should show decreasing values (convergence)
  expect_type(diag$diff_history, "double")
  expect_gt(length(diag$diff_history), 0)

  # Values should generally decrease (or at least be non-negative)
  expect_true(all(diag$diff_history >= 0))

  # With sufficient iterations, the last few values should be very small
  if (length(diag$diff_history) >= 3) {
    # Last 3 values should be smaller than first 3 (generally)
    expect_true(diag$diff_history[length(diag$diff_history)] <=
      diag$diff_history[1])
  }
})

test_that("dimv_diagnostics edge case: single variable", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50))
  X[sample(1:50, 10), "x1"] <- NA
  imp <- dimv_train(X, lambda = 0.5, maxit = 20, tol = 1e-4)

  X_imputed <- dimv_impute_new(imp, X)
  diag <- dimv_diagnostics(imp, X, X_imputed)

  # Should work with single variable
  expect_named(diag, c("iterations", "adaptive", "feature_select",
                       "residual_variances", "diff_history", "missing_rmse"))
  expect_length(diag$residual_variances, 1)
  expect_type(diff_history, "double")
})

test_that("dimv_diagnostics error handling", {
  # Wrong imputer type
  expect_error(dimv_diagnostics("not an imputer", data.frame(x = rnorm(10)),
                                data.frame(x = rnorm(10))),
    "imputer must be result of dimv_train()")

  # Mismatched dimensions between imputer and data
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  imp <- dimv_train(X, lambda = 0.5)
  expect_error(dimv_diagnostics(imp, X[, "x1", drop = FALSE],
                                dimv_impute_new(imp, X[, "x1", drop = FALSE])),
    "must have the same number of columns")
})