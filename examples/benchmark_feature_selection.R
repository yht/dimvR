#!/usr/bin/env Rscript

# Benchmark: Feature selection methods (lightweight)

suppressPackageStartupMessages(library(dimvR))

set.seed(42)

X <- as.data.frame(matrix(rnorm(500 * 30), ncol = 30))
colnames(X) <- paste0("V", seq_len(ncol(X)))
X$V2 <- 0.8 * X$V1 + rnorm(500, sd = 0.2)
X$V3 <- 0.6 * X$V1 + rnorm(500, sd = 0.3)

bench <- function(method, threshold = NULL) {
  t <- system.time({
    if (is.null(threshold)) {
      res <- select_features_adaptive(X, target_var = "V1", method = method)
    } else {
      res <- select_features_adaptive(X, target_var = "V1", method = method, threshold = threshold)
    }
  })
  list(
    method = method,
    elapsed = unname(t["elapsed"]),
    selected = length(res$selected_features)
  )
}

results <- rbind(
  as.data.frame(bench("adaptive")),
  as.data.frame(bench("fixed", threshold = 0.4)),
  as.data.frame(bench("mi")),
  as.data.frame(bench("hybrid"))
)

print(results)
