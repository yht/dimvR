# ===========================================
# cleanup.R
# Full cleanup utility for dimvExplainRpaper
# - removes experiment results, reports, caches, and future workers
# ===========================================

cleanup_project <- function(keep_reports = FALSE, verbose = TRUE, kill_futures = TRUE) {
  cat("🧹 Starting full cleanup of dimvExplainRpaper project...\n")
  
  root <- getwd()
  
  # stop parallel workers
  if (kill_futures) {
    if (requireNamespace("future", quietly = TRUE)) {
      tryCatch({
        future::plan(future::sequential)
        future::ClusterRegistry("stop")
      }, error = function(e) {})
    }
    if (requireNamespace("parallel", quietly = TRUE)) {
      try(parallel::stopCluster(parallel::makeCluster(0)), silent = TRUE)
    }
  }
  
  # target files to delete
  targets <- c(
    "dimv_results.RData",
    "dimv_results_pooled.RData",
    "dimv_full_report.html",
    "dimv_r_report.html",
    "Rplots.pdf"
  )
  
  if (!keep_reports) {
    pngs <- list.files(root, pattern = "\\.png$", full.names = TRUE)
    targets <- c(targets, pngs)
  }
  
  dirs_to_remove <- c(
    "src-i386", "src-x64", "build",
    "inst/doc", "tests", "tmp",
    "__pycache__", ".pytest_cache",
    ".Rproj.user", ".future"
  )
  
  deleted_files <- c()
  for (f in targets) {
    if (file.exists(f)) {
      file.remove(f)
      deleted_files <- c(deleted_files, f)
    }
  }
  
  for (d in dirs_to_remove) {
    if (dir.exists(d)) {
      unlink(d, recursive = TRUE, force = TRUE)
      deleted_files <- c(deleted_files, paste0(d, "/"))
    }
  }
  
  # clear temp directory used by R (session temporary files)
  temp_dirs <- list.dirs(tempdir(), recursive = FALSE, full.names = TRUE)
  if (length(temp_dirs)) {
    unlink(temp_dirs, recursive = TRUE, force = TRUE)
    deleted_files <- c(deleted_files, temp_dirs)
  }
  
  # clear cached workers from future
  if (requireNamespace("future", quietly = TRUE)) {
    future_cache <- Sys.getenv("R_FUTURE_CACHE_PATH", unset = NA)
    if (!is.na(future_cache) && dir.exists(future_cache)) {
      unlink(future_cache, recursive = TRUE, force = TRUE)
      deleted_files <- c(deleted_files, future_cache)
    }
  }
  
  # clean reticulate cache (Python venvs / .pyc)
  reticulate_dirs <- c("inst/python/__pycache__", "inst/python/.pytest_cache")
  for (rdir in reticulate_dirs) {
    if (dir.exists(rdir)) {
      unlink(rdir, recursive = TRUE, force = TRUE)
      deleted_files <- c(deleted_files, paste0(rdir, "/"))
    }
  }
  
  # remove .RData, .Rhistory
  cache_files <- c(".RData", ".Rhistory")
  for (cf in cache_files) {
    if (file.exists(cf)) {
      file.remove(cf)
      deleted_files <- c(deleted_files, cf)
    }
  }
  
  # remove any leftover compiled objects
  dll_files <- list.files("src", pattern = "\\.(o|so|dll)$", full.names = TRUE)
  if (length(dll_files)) {
    unlink(dll_files, force = TRUE)
    deleted_files <- c(deleted_files, dll_files)
  }
  
  # optional: devtools DLL cleaning
  if (requireNamespace("devtools", quietly = TRUE)) {
    try(devtools::clean_dll(), silent = TRUE)
  }
  
  # final report
  if (verbose) {
    if (length(deleted_files) == 0) {
      cat("✅ No temporary files found — project already clean.\n")
    } else {
      cat("✅ Cleanup complete. Deleted:\n")
      cat(paste0("   • ", deleted_files, collapse = "\n"), "\n")
    }
  }
  
  invisible(deleted_files)
}

# Auto-run when called directly via Rscript
if (sys.nframe() == 0L) {
  cleanup_project(keep_reports = FALSE)
}
