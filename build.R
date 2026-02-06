# Run this in R console to clean up and retry:

# 1) Remove lock/check artifacts (ignore errors if not exist)
try(unlink("dimvR.Rcheck/00LOCK-dimvR", recursive = TRUE, force = TRUE), silent = TRUE)
try(unlink("dimvR.Rcheck", recursive = TRUE, force = TRUE), silent = TRUE)

# 2) Clean and rebuild namespace/docs
try(file.remove("NAMESPACE"), silent = TRUE)
devtools::clean_dll()        # harmless even if no compiled code
devtools::document()

# 3) Run check with relaxed time/timestamp validations
withr::with_envvar(
  c(
    "_R_CHECK_TIMEZONE_" = "FALSE",
    "_R_CHECK_FUTURE_FILE_TIMESTAMPS_" = "FALSE",
    "_R_CHECK_SYSTEM_CLOCK_" = "FALSE"
  ),
  devtools::check(args = c("--ignore-vignettes", "--no-manual", "--as-cran"))
)
