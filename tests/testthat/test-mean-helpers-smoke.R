test_that("fit_mean_imputer computes means and handles non-data.frame input", {
  X <- data.frame(
    x1 = c(1, NA, 3),
    x2 = c(2, 4, NA)
  )

  imp <- dimvR:::fit_mean_imputer(X)
  expect_true(is.list(imp))
  expect_true(all(c("means", "cols") %in% names(imp)))
  expect_equal(unname(imp$means["x1"]), 2)
  expect_equal(unname(imp$means["x2"]), 3)

  expect_error(dimvR:::fit_mean_imputer(matrix(1:4, ncol = 2)), "data.frame")
})

test_that("apply_mean_imputer fills missing values using trained means", {
  X_train <- data.frame(
    x1 = c(1, 2, 3),
    x2 = c(4, 5, 6)
  )
  imp <- dimvR:::fit_mean_imputer(X_train)

  X_new <- data.frame(
    x1 = c(NA, 2),
    x2 = c(4, NA),
    x3 = c(10, 11)
  )

  X_imp <- dimvR:::apply_mean_imputer(imp, X_new)
  expect_equal(X_imp$x1[1], 2)
  expect_equal(X_imp$x2[2], 5)
  expect_equal(X_imp$x3, X_new$x3)

  expect_error(dimvR:::apply_mean_imputer(list(), X_new), "Invalid imputer")
})
