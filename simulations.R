source("functions.R")

normal_efficiency_sim <- function(N_pts, small_var = 1, estimator = "MM", initial_estimator = "S") {
    
    stopifnot(N_pts > 0, small_var > 0)
    
    if (estimator == "MM") {
        weight_fn <- function(x) bisquare_weight(x, k=4.68)
    } else if (estimator == "Huber") {
        initial_estimator <- "madn"
        weight_fn <- function(x) huber_weight(x)
    } else {
        stopifnot(1 == 0)
    }

    x_vals = seq(1, N_pts)
    
    small_noise <- rnorm(n = N_pts, mean = 0, sd = sqrt(small_var))
    y_vals = 1 * x_vals + small_noise

    bi_fit <- IRWLS_fit_simple(x_vals, y_vals, initial_est = initial_estimator, weight_fn = weight_fn)
    ols_fit <- lm(y_vals ~ x_vals)


    list(
        bisquare_slope = coef(bi_fit)[2],
        ols_slope      = coef(ols_fit)[2]
    )
}

