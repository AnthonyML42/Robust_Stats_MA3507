bisquare_rho <- function(x, k = 1.345) {
    ifelse(abs(x) <= k, 1 - (1 - (x/k)^2)^3, 1)
}


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
    stopifnot(length(x) == length(y), max_iter >= 0)

    X <- cbind(1, x)

    fit <- lm.fit(X, y)

    sigma_initial <- madn(residuals(fit))

    if (sigma_initial == 0) {
        med <- median(residuals(fit))
        if (med == 0) {
            print("ERROR")
            # maybe try eps as sigma_initial?
            stopifnot(1 == 0)
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


fast_s_estimator <- function(x, y, weight_fn = bisquare_weight, alpha = 0.01, p = 2, epsilon = 0.5) {
    # A Fast Algorithm for S-Regression Estimates 2006

    stopifnot(length(x) == length(y))

    N_subsamp_approx <- ceiling(-log(alpha) / ((1 - epsilon) ^ p)) # this uses an approximation of log on the denominator vs the paper

    for (i in 1:N_subsamp_approx) {
        subsample_idx <- sample(1:length(x), size=p)

        X <- cbind(1, x)

        fit <- lm.fit(X[subsample_idx, ], y[subsample_idx])

        resid <- y - X %*% fit$coefficients


        k <- 10 # apply k I-steps
        for (i in 1:k) {
            # 1. Compute the S-scale here and solve an optimisation problem with uniroot.
            # 2. Weighted least squares using the scale.
        }
        # 3. Compute the M-scale of all of them and keep the t best (store both their slope and scale estimates).
        # 4. Apply I-steps until the estimates all converge.
        # 5. Pick the one with the smallest scale and return it.
 
        stopifnot(1==0)
    }
}