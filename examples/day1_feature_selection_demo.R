#!/usr/bin/env Rscript

# Demo: Feature selection methods for dimvR

suppressPackageStartupMessages(library(dimvR))

set.seed(123)

# Simulated dataset
X <- data.frame(
  x1 = rnorm(200),
  x2 = rnorm(200),
  x3 = rnorm(200),
  x4 = rnorm(200),
  x5 = rnorm(200)
)

# Inject correlation
X$x2 <- 0.9 * X$x1 + rnorm(200, sd = 0.2)
X$x3 <- 0.6 * X$x1 + rnorm(200, sd = 0.4)

cat("=== Feature Selection Demo ===\n")
cat("Target: x1\n\n")

methods <- c("adaptive", "fixed", "mi", "hybrid")

for (m in methods) {
  if (m == "fixed") {
    res <- select_features_adaptive(X, target_var = "x1", method = m, threshold = 0.5)
  } else {
    res <- select_features_adaptive(X, target_var = "x1", method = m)
  }
  cat("[", m, "] selected:", paste(res$selected_names, collapse = ", "), "\n")
}

cat("\n=== High-dimensional demo ===\n")
X_hd <- as.data.frame(matrix(rnorm(100 * 50), ncol = 50))
colnames(X_hd) <- paste0("V", seq_len(ncol(X_hd)))
X_hd$V2 <- X_hd$V1 + rnorm(100, sd = 0.1)
X_hd$V3 <- 0.7 * X_hd$V1 + rnorm(100, sd = 0.3)

res_hd <- select_features_adaptive(X_hd, target_var = "V1", method = "hybrid", max_features = 5)
cat("Hybrid selected (HD):", paste(res_hd$selected_names, collapse = ", "), "\n")
