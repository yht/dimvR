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
  # Check multiple rows to increase chance of detecting variation
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
  # Test with m=1 (should be deterministic except for the noise, but with seed it should be reproducible)
  imputed_list_1 <- dimv_impute_multiple(imp, X, m = 1, seed = 456)
  imputed_list_1_again <- dimv_impute_multiple(imp, X, m = 1, seed = 456)
  expect_equal(imputed_list_1[[1]], imputed_list_1_again[[1]])
  # Test error handling
  expect_error(dimv_impute_multiple("not an imputer", X), "imputer must be result of dimv_train()")
})
