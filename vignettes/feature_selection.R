## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(dimvR)
set.seed(123)

X <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = rnorm(100),
  x4 = rnorm(100)
)

X$x2 <- 0.9 * X$x1 + rnorm(100, sd = 0.2)
X$x3 <- 0.6 * X$x1 + rnorm(100, sd = 0.3)

res_adaptive <- select_features_adaptive(X, target_var = "x1", method = "adaptive")
res_fixed <- select_features_adaptive(X, target_var = "x1", method = "fixed", threshold = 0.4)
res_mi <- select_features_adaptive(X, target_var = "x1", method = "mi")
res_hybrid <- select_features_adaptive(X, target_var = "x1", method = "hybrid")

res_adaptive$selected_names

## -----------------------------------------------------------------------------
imp <- dimv_train(X, lambda = 0.1, feature_select = TRUE)
print(imp)

