library(testthat)
library(dimvR)

test_that("dimv_impute_multiple returns correct structure and values", {
  set.seed(42)
  X <- data.frame(
    a = rnorm(60),
    b = rnorm(60),
    c = rnorm(60)
  )
  X[sample(1:nrow(X), 8), "a"] <- NA
  X[sample(1:nrow(X), 6), "b"] <- NA
  X[sample(1:nrow(X), 5), "c"] <- NA
  imp <- dimv_train(X, lambda = 0.2, maxit = 20, tol = 1e-4)
  imputed_list <- dimv_impute_multiple(imp, X, m = 5, seed = 123)
  expect_type(imputed_list, "list")
  expect_length(imputed_list, 5)
  for (i in seq_along(imputed_list)) {
    expect_s3_class(imputed_list[[i]], "data.frame")
    expect_equal(nrow(imputed_list[[i]]), nrow(X))
    expect_equal(ncol(imputed_list[[i]]), ncol(X))
    expect_equal(names(imputed_list[[i]]), names(X))
    expect_false(anyNA(imputed_list[[i]]))
  }
  # Check that the imputed values are different across imputations (due to added noise)
  diff_found <- FALSE
  for (j in seq_len(ncol(X))) {
    # Check first 5 rows
    col_vals <- sapply(imputed_list, function(df) df[[j]][1:5])
    if (length(unique(as.vector(col_vals))) > 1) {
      diff_found <- TRUE
      break
    }
  }
  expect_true(diff_found, info = "At least one column should show variation across imputations due to added noise")
})

test_that("dimv_impute_multiple works with various m values", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5)

  # Test m = 1 (deterministic with seed)
  imp1 <- dimv_impute_multiple(imp, X, m = 1, seed = 123)
  expect_length(imp1, 1)

  # Test m = 3
  imp3 <- dimv_impute_multiple(imp, X, m = 3, seed = 123)
  expect_length(imp3, 3)

  # Test m = 5
  imp5 <- dimv_impute_multiple(imp, X, m = 5, seed = 123)
  expect_length(imp5, 5)

  # Test m = 10
  imp10 <- dimv_impute_multiple(imp, X, m = 10, seed = 123)
  expect_length(imp10, 10)

  # All imputed datasets should have same dimensions
  for (lst in list(imp1, imp3, imp5, imp10)) {
    for (k in seq_along(lst)) {
      expect_equal(nrow(lst[[k]]), 50)
      expect_equal(ncol(lst[[k]]), 3)
      expect_equal(names(lst[[k]]), c("x1", "x2", "x3"))
    }
  }
})

test_that("dimv_impute_multiple is seed-reproducible", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5)

  # Run twice with same seed should give identical results
  imp_run1 <- dimv_impute_multiple(imp, X, m = 3, seed = 123)
  imp_run2 <- dimv_impute_multiple(imp, X, m = 3, seed = 123)

  for (i in seq_along(imp_run1)) {
    expect_equal(imp_run1[[i]], imp_run2[[i]],
      info = paste("Imputation", i, "should be identical with same seed")
    )
  }
})

test_that("dimv_impute_multiple with no missing values", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  imp <- dimv_train(X, lambda = 0.5)

  # No missing values - should still work
  imp_no_missing <- dimv_impute_multiple(imp, X, m = 3, seed = 123)
  expect_length(imp_no_missing, 3)
  for (i in seq_along(imp_no_missing)) {
    expect_false(anyNA(imp_no_missing[[i]]))
    expect_equal(nrow(imp_no_missing[[i]]), 50)
    expect_equal(names(imp_no_missing[[i]]), c("x1", "x2", "x3"))
  }
})

test_that("dimv_impute_multiple edge case: single column", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50))
  X[sample(1:nrow(X), 10), "x1"] <- NA
  imp <- dimv_train(X, lambda = 0.5)

  imp_single <- dimv_impute_multiple(imp, X, m = 3, seed = 123)
  expect_length(imp_single, 3)
  for (i in seq_along(imp_single)) {
    expect_equal(nrow(imp_single[[i]]), 50)
    expect_equal(names(imp_single[[i]]), "x1")
    expect_false(anyNA(imp_single[[i]]$x1))
  }
})

test_that("dimv_impute_multiple error handling", {
  # Wrong imputer type
  expect_error(dimv_impute_multiple("not an imputer", data.frame(x = rnorm(10))),
    "imputer must be result of dimv_train()")

  # m must be positive integer
  set.seed(123)
  X <- data.frame(x1 = rnorm(50))
  for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
  imp <- dimv_train(X, lambda = 0.5)

  # m = 0 should error or handle gracefully
  expect_error(dimv_impute_multiple(imp, X, m = 0),
    "m must be a positive integer")
})