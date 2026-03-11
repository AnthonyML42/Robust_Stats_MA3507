library(MASS)
source("functions.R")

# note we are currently not returning the intercept term!

normal_efficiency_sim <- function(N_pts, small_var = 1, estimator = "MM", initial_estimator = "S") {
    stopifnot(N_pts > 0, small_var > 0)

    if (estimator == "MM") {
        weight_fn <- function(x) bisquare_weight(x, k = 4.68)
    } else if (estimator == "Huber") {
        initial_estimator <- "madn"
        weight_fn <- function(x) huber_weight(x)
    } else {
        stopifnot(1 == 0)
    }

    # Generate evenly-spaced x values, then y=x+N(0, small_var) and fit OLS vs MM-estimator
    x_vals <- seq(1, N_pts)

    small_noise <- rnorm(n = N_pts, mean = 0, sd = sqrt(small_var))
    y_vals <- 1 * x_vals + small_noise

    bi_fit <- IRWLS_fit_simple(x_vals, y_vals, initial_est = initial_estimator, weight_fn = weight_fn)
    ols_fit <- lm(y_vals ~ x_vals)
    mass_fit <- rlm(y_vals ~ x_vals, method= "MM", psi = psi.bisquare)

    # return the slope [2], [1] would be the intercept
    list(
        bisquare_slope = coef(bi_fit)[2],
        ols_slope      = coef(ols_fit)[2],
        mass_slope     = coef(mass_fit)[2]
    )
}


gross_outlier_sim <- function(N_pts, small_var = 1, big_var = 100, contamination_proportion, estimator = "MM", initial_estimator = "S") {
    # TODO: How to let the caller decide the big noise distribution?
    # puts the job of making randomness on caller, current approach (randomness inside this function) more intuitive with furrr

    stopifnot(N_pts > 0, small_var > 0, contamination_proportion <= 1, contamination_proportion > 0)

    if (estimator == "MM") {
        weight_fn <- function(x) bisquare_weight(x, k = 4.68)
    } else if (estimator == "Huber") {
        initial_estimator <- "madn"
        weight_fn <- function(x) huber_weight(x)
    } else {
        stopifnot(1 == 0)
    }

    # this way with n_contaminates was chosen for easier file naming in run_simulations.R
    n_contaminates <- floor(contamination_proportion * N_pts)

    # Generate evenly spaced x values, then y+x+N(0, small_var)

    x_vals <- seq(1, N_pts)
    small_noise <- rnorm(n = N_pts, mean = 0, sd = sqrt(small_var))
    big_noise <- rnorm(n = n_contaminates, mean = 0, sd = sqrt(big_var))

    y_vals <- 1 * x_vals + small_noise

    # Corrupt a proportion by adding noise to them and then fit OLS vs MM-estimator
    contamination_idx <- sample(1:N_pts, size = n_contaminates)

    y_vals[contamination_idx] <- y_vals[contamination_idx] + big_noise

    bi_fit <- IRWLS_fit_simple(x_vals, y_vals, initial_est = initial_estimator, weight_fn = weight_fn)
    ols_fit <- lm(y_vals ~ x_vals)
    mass_fit <- rlm(y_vals ~ x_vals, method= "MM", psi = psi.bisquare)

    list(
        bisquare_slope = coef(bi_fit)[2],
        ols_slope      = coef(ols_fit)[2],
        mass_slope     = coef(mass_fit)[2]
    )
}
