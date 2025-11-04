library(dimvExplainRpaper)
library(future)
library(magrittr)
plan(multisession)

# prepare dataset (California as before)
caldata <- data.table::fread("https://raw.githubusercontent.com/ageron/handson-ml/master/datasets/housing/housing.csv")
cal <- caldata %>% dplyr::select(-ocean_proximity) %>% head(100)
X_cal <- cal %>% dplyr::select(-median_house_value) %>% as.data.frame()
X_cal[] <- lapply(X_cal, as.numeric)
y_cal <- cal$median_house_value

# run with pooling
res <- run_full_pipeline(X_cal, y_cal,
                         task = "regression",
                         missing_rates = c(0.2, 0.4),
                         mechanisms = c("MCAR", "MNAR"),
                         imputers = c("dimv_R", "mice", "mean"),
                         nsim = 500,
                         repeats = 1,
                         workers = 4,
                         m_imp = 5,
                         seed = 2025)

save(res, file = "dimv_results_pooled.RData")
generate_report(res, out_file = "dimv_full_report_pooled.html", dataset_name = "California Housing")
browseURL("dimv_full_report_pooled.html")
