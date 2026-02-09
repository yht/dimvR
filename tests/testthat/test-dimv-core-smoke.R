test_that("dimv_train returns a valid imputer object", {
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

  expect_s3_class(imp, "dimv_imputer")
  expect_equal(length(imp$models), ncol(X))
  expect_equal(length(imp$col_means), ncol(X))
  expect_true(is.numeric(imp$resid_var))
})

test_that("dimv_impute_new imputes all missing values and preserves schema", {
  set.seed(99)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  X[sample(1:nrow(X), 7), "x1"] <- NA
  X[sample(1:nrow(X), 7), "x2"] <- NA
  X[sample(1:nrow(X), 7), "x3"] <- NA

  imp <- dimv_train(X, lambda = 0.1, maxit = 20, tol = 1e-4)
  X_imp <- dimv_impute_new(imp, X, maxit = 20, tol = 1e-4)

  expect_equal(names(X_imp), names(X))
  expect_equal(nrow(X_imp), nrow(X))
  expect_false(anyNA(X_imp))
})
