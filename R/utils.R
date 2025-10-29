# ===========================================
# File: R/utils.R
# ===========================================

#' Calculate Mean Squared Error
#' 
#' @param a Numeric vector or convertible to numeric
#' @param b Numeric vector or convertible to numeric
#' @return Mean squared error between a and b
#' @keywords internal
mse <- function(a, b) {
  mean((as.numeric(a) - as.numeric(b))^2, na.rm = TRUE)
}

# Internal cache environment for imputer objects
.imputer_cache <- new.env(parent = emptyenv())

#' Generate Cache Key
#' 
#' @param tag Character string tag
#' @param X Data object to hash
#' @return Cache key string
#' @keywords internal
cache_key <- function(tag, X) {
  paste0(tag, "_", digest::digest(X))
}

#' Retrieve Cached Object
#' 
#' @param tag Character string tag
#' @param X Data object (for key generation)
#' @return Cached object or NULL if not found
#' @keywords internal
get_cached <- function(tag, X) {
  key <- cache_key(tag, X)
  if (exists(key, envir = .imputer_cache)) {
    get(key, envir = .imputer_cache)
  } else {
    NULL
  }
}

#' Store Object in Cache
#' 
#' @param tag Character string tag
#' @param X Data object (for key generation)
#' @param obj Object to cache
#' @return NULL (invisibly)
#' @keywords internal
set_cached <- function(tag, X, obj) {
  key <- cache_key(tag, X)
  assign(key, obj, envir = .imputer_cache)
  invisible(NULL)
}