fit_mean_imputer <- function(X_train) {
  list(means = sapply(X_train, function(x) mean(x, na.rm = TRUE)))
}

apply_mean_imputer <- function(imputer, X) {
  as.data.frame(mapply(function(col, m) { col[is.na(col)] <- m; col }, X, imputer$means, SIMPLIFY=FALSE))
}
