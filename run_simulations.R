source("simulations.R")

library(furrr)
library(purrr)
library(dplyr)
library(readr)
library(glue)

normal_eff_outpath <- "data/normal_efficiency.csv"


plan(multisession, workers = parallelly::availableCores() - 1)

if (!file.exists(normal_eff_outpath)) {
    cat("Starting normal efficiency simulation...\n")
    future_map(
        seq_len(10000),
        \(rep_i) map(c(50, 200, 500), \(n) {
            sim <- normal_efficiency_sim(N_pts = n, small_var = 1)
            tibble(rep = rep_i, N_pts = n, slope_mm = sim$bisquare_slope, slope_ols = sim$ols_slope, intercept_mm = sim$bisquare_int, intercept_ols = sim$ols_int)
        }) %>% list_rbind(),
        .options = furrr_options(seed = 0),
        .progress = TRUE
    ) %>%
        list_rbind() %>%
        write_csv(normal_eff_outpath)
} else {
    cat("Normal efficiency file found. Skipping...\n")
}

contamination_pct <- 40
gross_outlier_outpath <- glue("data/gross_outlier_{contamination_pct}.csv")

if (!file.exists(gross_outlier_outpath)) {
    cat("Starting gross outlier simulation...\n")
    future_map(
        seq_len(10000),
        \(rep_i) map(c(50, 200, 500), \(n) {
            sim <- gross_outlier_sim(N_pts = n, small_var = 1, big_var = 100, contamination_proportion = contamination_pct / 100)
            tibble(rep = rep_i, N_pts = n, slope_mm = sim$bisquare_slope, slope_ols = sim$ols_slope, intercept_mm = sim$bisquare_int, intercept_ols = sim$ols_int)
        }) %>% list_rbind(),
        .options = furrr_options(seed = 0),
        .progress = TRUE
    ) %>%
        list_rbind() %>%
        write_csv(gross_outlier_outpath)
} else {
    cat("Gross outlier file found. Skipping...\n")
}


breakdown_outpath <- "data/breakdown.csv"
if (!file.exists(breakdown_outpath)) {
    cat("Starting breakdown simulation...\n")
    future_map(
        seq_len(1500),
        \(rep_i) map(c(0.2, 0.45, 0.49, 0.51, 0.55), \(c_p) {
            sim <- gross_outlier_sim(N_pts = 200, small_var = 1, big_mean = 10000, big_var = 1, contamination_proportion = c_p)
            tibble(rep = rep_i, N_pts = 200, slope_mm = sim$bisquare_slope, slope_ols = sim$ols_slope, intercept_mm = sim$bisquare_int, intercept_ols = sim$ols_int, cp = c_p)
        }) %>% list_rbind(),
        .options = furrr_options(seed = 0),
        .progress = TRUE
    ) %>%
        list_rbind() %>%
        write_csv(breakdown_outpath)
} else {
    cat("Breakdown file found. Skipping...\n")
}
plan(sequential)
