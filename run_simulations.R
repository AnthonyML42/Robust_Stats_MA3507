source("simulations.R")

library(furrr)
library(purrr)
library(dplyr)
library(readr)

plan(multisession, workers = 8)

future_map(
    seq_len(10000),
    \(rep_i) map(c(50, 200, 500), \(n) {
        sim <- normal_efficiency_sim(N_pts = n, small_var = 1)
        tibble(rep = rep_i, N_pts = n, slope_mm = sim$bisquare_slope, slope_ols = sim$ols_slope)
    }) %>% list_rbind(),
    .options = furrr_options(seed = 0),
    .progress = TRUE
) %>%
    list_rbind() %>%
    write_csv("data/normal_efficiency.csv")


plan(sequential)
