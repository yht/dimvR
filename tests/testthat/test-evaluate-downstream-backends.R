library(testthat)
library(dimvR)

make_backend_data <- function(n = 80) {
  set.seed(101)
  data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
}

test_that("evaluate_downstream uses identical formula columns for train and test", {
  dat <- make_backend_data()
  fit <- evaluate_downstream(dat, y ~ x1 + x2, model_method = "lm",
                             test_idx = 1:20)

  expect_s3_class(fit$model, "lm")
  expect_length(fit$predictions, 20)
  expect_true(is.finite(fit$mse))
  expect_equal(fit$method, "lm")
})

test_that("optional downstream backends return predictions", {
  dat <- make_backend_data()
  methods <- c("ranger", "randomForest", "glmnet")

  for (method in methods) {
    skip_if_not_installed(method)
    fit <- evaluate_downstream(dat, y ~ x1 + x2, model_method = method,
                               test_idx = 1:20)
    expect_length(fit$predictions, 20)
    expect_true(is.finite(fit$mse))
    expect_equal(fit$method, method)
  }
})

test_that("missing optional backend has an actionable error", {
  skip_if(requireNamespace("ranger", quietly = TRUE))
  dat <- make_backend_data()
  expect_error(
    evaluate_downstream(dat, y ~ x1 + x2, model_method = "ranger"),
    "Package 'ranger' needed"
  )
})

test_that("evaluator supports matrix input, default split, custom models, and glmnet s", {
  dat <- make_backend_data()
  mat <- as.matrix(dat)
  fit <- evaluate_downstream(mat, "y ~ x1 + x2", model_method = "lm")
  expect_length(fit$predictions, 16)

  prefit <- lm(y ~ x1 + x2, data = dat)
  custom <- evaluate_downstream(dat, y ~ x1 + x2, model_method = "custom",
                                model = prefit, test_idx = 1:10)
  expect_equal(custom$method, "custom")
  expect_length(custom$predictions, 10)

  skip_if_not_installed("glmnet")
  penalized <- evaluate_downstream(dat, y ~ x1 + x2, model_method = "glmnet",
                                   test_idx = 1:10, s = 0.01)
  expect_length(penalized$predictions, 10)
})
