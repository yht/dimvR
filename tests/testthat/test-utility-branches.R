library(testthat)
library(dimvR)

test_that("diagnostic helpers score features and select lambda levels", {
  set.seed(11)
  z <- rnorm(60)
  X <- data.frame(a = z, b = z + rnorm(60, sd = .05), c = rnorm(60))
  scores <- feature_select_score(X)
  expect_equal(scores$feature, names(X))
  expect_true(all(is.finite(scores$score)))
  expect_equal(adaptive_lambda(data.frame(a = z, b = z + rnorm(60, sd = .01))), 5)
  expect_equal(adaptive_lambda(data.frame(a = z, b = z + rnorm(60, sd = .2))), 5)
  expect_equal(adaptive_lambda(data.frame(a = z, b = z + rnorm(60, sd = .8))), 2)
  expect_equal(adaptive_lambda(data.frame(a = z, b = rnorm(60))), .5)
})

test_that("diagnostic summary, print, and plot helpers work", {
  set.seed(12)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  X$a[1:4] <- NA
  imp <- dimv_train(X, maxit = 5)
  d <- dimv_convergence_diag(imp)
  expect_equal(d$iterations, imp$iters)
  expect_equal(d$lambda, imp$lambda)
  expect_output(print_dimv_diag(imp), "DIMV Diagnostics")
  expect_error(dimv_convergence_diag(NULL), "Requires a dimv_imputer")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(plot_feature_selection(feature_select_score(X)), "ggplot")
  }
})

test_that("experiment utility functions cover pooling and missingness branches", {
  expect_equal(dimvR:::within_var_mse(1, 1), 0)
  expect_gt(dimvR:::within_var_mse(1:3, c(0, 2, 4)), 0)
  expect_equal(dimvR:::within_var_shap_mse(matrix(1), matrix(1)), 0)
  expect_gt(dimvR:::within_var_shap_mse(matrix(1:4, 2), matrix(0, 2, 2)), 0)
  expect_error(dimvR:::within_var_shap_mse(matrix(1), matrix(1, 1, 2)), "same shape")

  one <- dimvR:::rubin_pool(1, 0.25)
  many <- dimvR:::rubin_pool(c(1, 2), c(.1, .2))
  expect_equal(one$m, 1)
  expect_equal(one$Q_bar, 1)
  expect_true(many$between > 0)
  expect_true(many$se > 0)

  set.seed(13)
  X <- data.frame(a = rnorm(40), b = rnorm(40), c = rnorm(40))
  for (mechanism in c("MCAR", "MAR", "MNAR")) {
    miss <- dimvR:::simulate_missing(X, rate = .3, mechanism = mechanism)
    expect_equal(dim(miss), dim(X))
    expect_true(anyNA(miss))
  }
  expect_error(dimvR:::simulate_missing(X, mechanism = "bad"), "Unknown mechanism")
  expect_error(dimvR:::.require_pkg("package_that_does_not_exist"), "required but not installed")
})

test_that("mean imputer validates inputs and aligns known columns", {
  expect_error(dimvR:::fit_mean_imputer(list(x = 1:3)), "data.frame")
  X <- data.frame(a = c(1, NA, 3), b = c("x", "y", "z"))
  expect_warning(fit <- dimvR:::fit_mean_imputer(X), "Non-numeric")
  expect_equal(unname(fit$means["a"]), 2)
  out <- dimvR:::apply_mean_imputer(fit, data.frame(b = NA, a = NA, extra = NA))
  expect_equal(out$a, 2)
  expect_true(is.na(out$b))
  expect_error(dimvR:::apply_mean_imputer(list(), X), "Invalid imputer")
})

test_that("feature-selection fallback and empty-candidate branches are covered", {
  one <- data.frame(x = rnorm(20))
  empty <- select_features_adaptive(one, target_var = 1, method = "adaptive")
  expect_length(empty$selected_features, 0)
  expect_error(select_features_adaptive(data.frame(x = one$x, y = rnorm(20)),
                                        target_var = 1, method = "fixed"),
               "threshold must be provided")
  expect_error(select_features_adaptive(one, target_var = 0), "valid column")

  scores <- c(a = .1, b = .2, c = .3)
  expect_length(dimvR:::select_adaptive_threshold(scores, min_features = 2)$selected, 2)
  expect_length(dimvR:::select_fixed_threshold(scores, threshold = .99, min_features = 2)$selected, 2)
  expect_length(dimvR:::select_adaptive_threshold(numeric(0))$selected, 0)
  expect_length(dimvR:::select_fixed_threshold(numeric(0))$selected, 0)

  X <- data.frame(a = 1:20, b = c(rep(NA, 18), 1, 2), c = rnorm(20))
  mi <- select_features_adaptive(X, target_var = 1, method = "mi", min_features = 1)
  hy <- select_features_adaptive(X, target_var = 1, method = "hybrid", min_features = 1)
  expect_true(is.list(mi$method_info))
  expect_true(is.list(hy$method_info))

  empty_scores <- numeric(0)
  expect_length(dimvR:::select_mutual_information(one, 1)$selected, 0)
  expect_length(dimvR:::select_hybrid(one, 1)$selected, 0)
  expect_length(dimvR:::select_fixed_threshold(empty_scores)$selected, 0)
})

test_that("DIMV optional feature selection and adaptive regularization paths work", {
  set.seed(15)
  z <- rnorm(40)
  X <- data.frame(a = z, b = z + rnorm(40, .1), c = rnorm(40))
  X$b[1:5] <- NA
  fit <- dimv_train(X, lambda = .2, maxit = 3, feature_select = TRUE,
                    adaptive = TRUE, verbose = TRUE)
  expect_s3_class(fit, "dimv_imputer")
  expect_true(all(vapply(fit$models, function(x) length(x$features) > 0, logical(1))))
  expect_gt(dimvR:::adaptive_lambda_fn(matrix(c(1, 1, 2, 2), 2)), 0)
  expect_length(dimvR:::feature_select_fn(as.matrix(X), 1, max_features = 2), 2)

  complete <- dimv_train(X[, c("a", "c")], maxit = 2, verbose = TRUE)
  expect_output(print(complete), "R-native DIMV imputer")
  expect_output(dimv_impute_new(complete, X[, c("a", "c")], verbose = TRUE),
                "no missing")
})

test_that("multiple imputation validates m and preserves observed values", {
  set.seed(16)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  X$a[1:4] <- NA
  imp <- dimv_train(X, maxit = 3)
  expect_error(dimv_impute_multiple(imp, X, m = 1.5), "positive integer")
  expect_error(dimv_impute_multiple("bad", X), "imputer must be")
  draws <- dimv_impute_multiple(imp, X, m = 2, seed = 16)
  expect_length(draws, 2)
  expect_equal(draws[[1]]$b, X$b)
})

test_that("MICE backend trains and applies train-only imputations", {
  skip_if_not_installed("mice")
  set.seed(14)
  train <- data.frame(a = c(1, NA, 3, 4, 2, 6, 7, 5, 9, 8, 10, 11),
                      b = c(6, 2, NA, 3, 5, 1, 8, 11, 4, 10, 7, 9))
  fit <- dimvR:::mice_backend_train(train, m = 2, maxit = 1, seed = 14)
  expect_length(fit$imputations, 2)
  test <- data.frame(a = c(NA, 8), b = c(2, NA))
  applied <- dimvR:::mice_backend_apply(fit, test, k = 1)
  expect_equal(dim(applied), dim(test))
  expect_false(anyNA(applied))
  expect_error(dimvR:::mice_backend_apply(fit, test, k = 3), "between 1 and")
})
