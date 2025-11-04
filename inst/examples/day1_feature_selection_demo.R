# ===========================================
# Phase 1, Day 1 Demo: Feature Selection
# ===========================================

library(dimvR)

# Set seed for reproducibility
set.seed(2024)

cat("============================================\n")
cat("DIMV Feature Selection Demo\n")
cat("Phase 1 - Day 1 Implementation\n")
cat("============================================\n\n")

# ========== Example 1: Basic Adaptive Selection ==========
cat("Example 1: Adaptive Feature Selection\n")
cat("--------------------------------------\n")

# Create synthetic dataset with known correlation structure
n <- 100
X <- data.frame(
  target = rnorm(n),
  high_cor = rnorm(n),
  medium_cor = rnorm(n),
  low_cor = rnorm(n),
  noise1 = rnorm(n),
  noise2 = rnorm(n)
)

# Create correlations
X$high_cor <- 0.9 * X$target + rnorm(n, sd = 0.3)
X$medium_cor <- 0.5 * X$target + rnorm(n, sd = 0.5)
X$low_cor <- 0.2 * X$target + rnorm(n, sd = 0.8)

cat("Dataset created with known correlation structure:\n")
cat("  - high_cor: r ~ 0.9\n")
cat("  - medium_cor: r ~ 0.5\n")
cat("  - low_cor: r ~ 0.2\n")
cat("  - noise1, noise2: r ~ 0\n\n")

# Test adaptive selection
result_adaptive <- select_features_adaptive(
  X, 
  target_var = "target",
  method = "adaptive",
  min_features = 2,
  verbose = TRUE
)

cat("\nSelected features (adaptive):\n")
print(colnames(X)[result_adaptive$selected_features])

cat("\nAll correlation scores:\n")
print(round(result_adaptive$scores, 3))


# ========== Example 2: Compare All Methods ==========
cat("\n\n============================================\n")
cat("Example 2: Method Comparison\n")
cat("--------------------------------------\n")

methods <- c("adaptive", "fixed", "mi", "hybrid")
results <- list()

for (method in methods) {
  cat("\nTesting method:", method, "\n")
  
  if (method == "fixed") {
    results[[method]] <- select_features_adaptive(
      X, 
      target_var = 1,
      method = method,
      threshold = 0.4,
      min_features = 2,
      verbose = FALSE
    )
  } else {
    results[[method]] <- select_features_adaptive(
      X, 
      target_var = 1,
      method = method,
      min_features = 2,
      verbose = FALSE
    )
  }
  
  cat("  Selected", length(results[[method]]$selected_features), "features:",
      paste(colnames(X)[results[[method]]$selected_features], collapse = ", "), "\n")
}


# ========== Example 3: With Missing Data ==========
cat("\n\n============================================\n")
cat("Example 3: Selection with Missing Data\n")
cat("--------------------------------------\n")

# Create copy with missing values
X_missing <- X
X_missing$high_cor[sample(n, 10)] <- NA
X_missing$medium_cor[sample(n, 10)] <- NA

cat("Introduced 10% missing data in two variables\n\n")

result_missing <- select_features_adaptive(
  X_missing,
  target_var = 1,
  method = "adaptive",
  verbose = TRUE
)

cat("\nSelected features with missing data:\n")
print(colnames(X_missing)[result_missing$selected_features])


# ========== Example 4: High-Dimensional Case ==========
cat("\n\n============================================\n")
cat("Example 4: High-Dimensional Selection\n")
cat("--------------------------------------\n")

# Create dataset with many variables
n <- 200
p <- 50
X_large <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(X_large) <- paste0("V", 1:p)

# Make V2-V5 correlated with V1
for (j in 2:5) {
  X_large[, j] <- 0.7 * X_large[, 1] + rnorm(n, sd = 0.5)
}

cat("Created dataset with", p, "variables\n")
cat("Variables V2-V5 are correlated with V1\n\n")

result_large <- select_features_adaptive(
  X_large,
  target_var = 1,
  method = "adaptive",
  max_features = 10,
  verbose = TRUE
)

cat("\nTop selected features:\n")
top_features <- colnames(X_large)[result_large$selected_features]
cat(paste(top_features, collapse = ", "), "\n")

cat("\nTop 10 correlation scores:\n")
top_scores <- sort(result_large$scores, decreasing = TRUE)[1:10]
print(round(top_scores, 3))


# ========== Example 5: Integration with DIMV ==========
cat("\n\n============================================\n")
cat("Example 5: Integration with DIMV Imputation\n")
cat("--------------------------------------\n")

# Use the original dataset with missing data
X_test <- X
X_test$target[sample(n, 15)] <- NA
X_test$high_cor[sample(n, 15)] <- NA

cat("Dataset with missing values:\n")
cat("  Missing in target:", sum(is.na(X_test$target)), "\n")
cat("  Missing in high_cor:", sum(is.na(X_test$high_cor)), "\n\n")

# Note: Future integration point
cat("Future integration:\n")
cat("  - DIMV will use adaptive feature selection internally\n")
cat("  - Each variable will be imputed using only relevant predictors\n")
cat("  - This reduces noise and improves convergence\n\n")

# For now, demonstrate what features would be used
for (var_idx in which(sapply(X_test, function(x) any(is.na(x))))) {
  var_name <- colnames(X_test)[var_idx]
  sel <- select_features_adaptive(X_test, target_var = var_idx, 
                                  method = "adaptive", verbose = FALSE)
  
  cat("For imputing", var_name, ":\n")
  cat("  Would use:", paste(colnames(X_test)[sel$selected_features], 
                           collapse = ", "), "\n")
}


# ========== Summary ==========
cat("\n\n============================================\n")
cat("Summary of Day 1 Implementation\n")
cat("============================================\n")
cat("\n✓ Implemented 4 feature selection methods:\n")
cat("  1. Adaptive (median-based threshold)\n")
cat("  2. Fixed threshold\n")
cat("  3. Mutual Information\n")
cat("  4. Hybrid (correlation + MI)\n")
cat("\n✓ Features:\n")
cat("  - Handles missing data automatically\n")
cat("  - Respects min/max feature constraints\n")
cat("  - Works with high-dimensional data\n")
cat("  - Ready for DIMV integration\n")
cat("\n✓ Next steps (Days 2-3):\n")
cat("  - Add comprehensive unit tests\n")
cat("  - Create benchmarking suite\n")
cat("  - Write documentation and vignette\n")
cat("  - Integrate with main DIMV algorithm\n")
cat("\n============================================\n\n")
