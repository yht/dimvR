# ===========================================
# File: R/utils.R
# ===========================================

mse <- function(a,b) mean((as.numeric(a)-as.numeric(b))^2, na.rm=TRUE)

if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
library(digest)

.imputer_cache <- new.env(parent = emptyenv())

cache_key <- function(tag, X) paste0(tag, "_", digest::digest(X))

get_cached <- function(tag, X) {
  key <- cache_key(tag, X)
  if (exists(key, envir = .imputer_cache)) get(key, envir = .imputer_cache) else NULL
}

set_cached <- function(tag, X, obj) {
  key <- cache_key(tag, X)
  assign(key, obj, envir = .imputer_cache)
}
