huber_psi <- function(x, k = 1.345) {
    # constant for 95% normal efficiency of the estimator

    return (ifelse(abs(x) <= k, x, k * sign(x)))
}

bisquare_psi <- function(x, k = 4.68) {
    ifelse(abs(x) <= k, (6 / k^2) * x * (1 - (x/k)^2)^2, 0)
}

bisquare_weight <- function(x, k = 4.68) {
    ifelse(abs(x) <= k, (1 - (x/k)^2)^2, 0)
}

huber_weight <- function(x, k = 1.345) {
    ifelse(abs(x) <= k, 1, k / abs(x))
}

madn <- function(x) {
    c = 0.675 # Phi^-1(0.75) normal dist

    med = median(x, na.rm = TRUE)
    return(1/c * median(abs(x - med), na.rm = TRUE))
}

IRWLS_fit_simple <- function(x, y, weight_fn, max_iter = 1000, eps = 1e-6) {
    # The Iteratively Reweighted Least Squares algorithm is used for fitting least squares with "robust" loss functions.
    # Simple because the dispersion is not estimated from an S-estimator, and the initial residuals are done with OLS. 
    # This can cause the estimate to converge to a local minima for nonconvex psi functions.

    # using lm naively is much slower vs lm.fit and the design matrix X.
    stopifnot(length(x) == length(y), max_iter <= 0)

    X <- cbind(1, x)

    fit <- lm.fit(X, y)

    sigma_initial <- madn(residuals(fit))

    if (sigma_initial == 0) {
        med <- median(residuals(fit))
        if (med == 0) {
            print("ERROR")
            # maybe try eps as sigma_initial?
            stopifnot(1==0)
        }
    }

    for (k in 1:max_iter) {
        resid <- residuals(fit)
        candidate_weights <- weight_fn(resid / sigma_initial)

        candidate_fit <- lm.wfit(X, y, candidate_weights) 

        if (max(abs(residuals(fit) - residuals(candidate_fit))) < eps) {
            break
        } else {
            fit <- candidate_fit
        }


    }

    return(fit)

}


fast_s_estimator <- function(x, y, alpha = 0.01, p = 2, epsilon = 0.5) {
    # from A Fast Procedure for Outlier Diagnostics in Large Regression Problems 1999
    N_subsamp_approx <- ceil(-log(alpha) / ((1 - epsilon) ^ p))
}