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

IRWLS_fit_simple <- function(x, y, weight_fn, initial_est = "madn", max_iter = 1000, eps = 1e-6) {
    # The Iteratively Reweighted Least Squares algorithm is used for fitting least squares with "robust" loss functions.
    # Simple because the dispersion is not estimated from an S-estimator, and the initial residuals are done with OLS. 
    # This can cause the estimate to converge to a local minima for nonconvex psi functions.

    # using lm naively is much slower vs lm.fit and the design matrix X.
    stopifnot(length(x) == length(y), max_iter >= 0)

    X <- cbind(1, x)


    if (initial_est == "madn") {
        fit <- lm.fit(X, y)
        sigma_initial <- madn(residuals(fit))
    } else if (initial_est == "S") {
        fast_s <- fast_s_estimator(x, y)
        # originally fast_s returned the beta vector and scale, but I want to use residuals(...) so it needs to return a fit
        sigma_initial <- fast_s$scale
        fit <- fast_s$fit
    } else {
        print("Initial estimate not implemented")
        stopifnot(1 == 0)
    }

    
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

s_scale_iteration <- function(r, start_val = madn(r), rho_fn, max_iter = 1000, eps = 1e-6) {
    # Estimation of scale according to 2.8.2 Maronna
    sigma_est <- start_val
    b <- 0.5
    for (i in 1:max_iter) {

        rho_vec <- as.vector(rho_fn(r / sigma_est))
        c = 1 / (length(r) * b)
        sigma_new <- sigma_est * sqrt(mean(rho_vec) / b)

        if (abs(sigma_new/sigma_est - 1) < eps) {
            sigma_est <- sigma_new
            break
        } else {
            sigma_est <- sigma_new
        }
    }
    return(sigma_est)
}


fast_s_estimator <- function(x, y, alpha = 0.01, p = 2, max_iter = 1500, eps = 1e-6) {
    # A Fast Algorithm for S-Regression Estimates 2006
    # Because we're focusing on the p=2 case, the line (between two points) is a perfect fit
    # This leads to some practical numerical problems, so we can increase the size of the subsample to 3 instead of 2 elements - and this doesn't affect the quality of the resulting estimator

    stopifnot(length(x) == length(y))

    # take the maximum breakdown epsilon=0.5 here, not to be confused with the convergence parameter "eps"
    N_subsample_approx <- ceiling(-log(alpha) / ((1 - 0.5) ^ p)) # this uses an approximation of log on the denominator vs the paper

    rho_fn <- function(x) bisquare_rho(x, k=1.547)
    weight_fn <- function(x) bisquare_weight(x, k=1.547)

    candidate_prelim <- vector("list", N_subsample_approx)
    X <- cbind(1, x)

    for (i in 1:N_subsample_approx) {
        subsample_idx <- sample(1:length(x), size=(p+1))

        # fit a line through the subsampled points and then take the residuals with this line on the whole data

        fit <- lm.fit(X[subsample_idx, ], y[subsample_idx])

        resid <- y - X %*% fit$coefficients
 
        sigma <- s_scale_iteration(resid, rho_fn=rho_fn)

        # Perform Improvement-steps on the candidates
        # the following is very similar to IRWLS, hence the name I-step. Maybe the IRWLS function can be adapted to be used here
        # IRWLS_fit_simple does not currently jointly estimate scale - which is needed in this case

        k <- 2 

        for (j in 1:k) {
            weights <- as.vector(weight_fn(resid / sigma))
            fit <- lm.wfit(X, y, weights)
            resid <- residuals(fit)
            sigma <- s_scale_iteration(resid, sigma, rho_fn=rho_fn)

        }

        candidate_prelim[[i]] <- list(beta = fit$coefficients, scale = sigma)
    }
    # Sort and only keep the "good" ones
    # How many to keep?

    N_keep <- min(10, N_subsample_approx) # this for now
    results <- vector("list", N_keep)

    sorted_scale_idx <- sapply(candidate_prelim, `[[`, "scale")
    candidate_prelim <- candidate_prelim[order(sorted_scale_idx)[1:N_keep]]


    # Now I-steps until convergence
    for (c in seq_along(candidate_prelim)) { # R doesn't have .enumerate()!

        beta <- as.vector(candidate_prelim[[c]]$beta)
        sigma <- candidate_prelim[[c]]$scale

        resid <- y - X %*% beta

        # I-steps until convergence
        for (i in 1:max_iter) {
            weights <- as.vector(weight_fn(resid / sigma))
            fit <- lm.wfit(X, y, weights)
            resid <- residuals(fit)
            sigma_new <- s_scale_iteration(resid, sigma, rho_fn=rho_fn)

            if (abs((sigma_new / sigma) - 1) < eps) {
                sigma <- sigma_new
                break
            } else {
                sigma <- sigma_new
            }

        }
        results[[c]] <- list(fit = fit, scale = sigma)
    }

    min_scale_est <- results[[which.min(sapply(results, `[[`, "scale"))]]

    return(list(fit = min_scale_est$fit, scale = min_scale_est$scale))

}


    #k <- 10 # apply k I-steps
    #for (i in 1:k) {
        # 1. Compute the S-scale here and used fixed point iteration - no uniroot needed!
        # 2. Weighted least squares using the scale.
    #}
    # 3. Compute the M-scale of all of them and keep the t best (store both their slope and scale estimates).
    # 4. Apply I-steps until the estimates all converge.
    # 5. Pick the one with the smallest scale and return it.
