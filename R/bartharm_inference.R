# This function runs the BARTharm algorithm for Bayesian Additive Regression Trees (BART) harmonization.
# It performs Gibbs sampling to estimate harmonization parameters.
#
# Arguments:
# - num_iter: Number of MCMC iterations
# - thinning_interval: Interval for thinning MCMC draws
# - X_iqm_matrix: Normalised matrix of IQM covariates
# - X_bio_matrix: Normalised matrix of biological covariates
# - Y: Outcome variable
# - hypers_mu, hypers_tau: Hyperparameter settings for mu and tau forests
# - opts_mu, opts_tau: Options for mu and tau forests
# - var_scaling: Logical. If TRUE, indicates that variance harmonization will be performed
# - site_labels: Vector of site/scanner IDs corresponding to each subject
# - alpha0, beta0: Hyperparameters for inverse-gamma prior on site-specific variances

bartharm_inference <- function(num_iter,
                                    thinning_interval,
                                    X_iqm_matrix,
                                    X_bio_matrix,
                                    Y,
                                    hypers_mu,
                                    hypers_tau,
                                    opts_mu,
                                    opts_tau,
                                    var_scaling = TRUE,
                                    site_labels = NULL,
                                    alpha0 = 0.01,
                                    beta0  = 0.01) {
  # Minimal-change Gibbs sampler. DOES NOT fix a reference site or renormalise deltas.
  n <- length(Y)
  num_saved_iters <- floor(num_iter / thinning_interval)
  if (num_saved_iters < 1) stop("thinning_interval too large or num_iter too small")
  
  # storage
  mu_out  <- matrix(NA, nrow = num_saved_iters, ncol = n)
  tau_out <- matrix(NA, nrow = num_saved_iters, ncol = n)
  sigma_out <- numeric(num_saved_iters)
  
  if (var_scaling) {
    if (is.null(site_labels)) stop("site_labels must be provided when var_scaling = TRUE")
    site_list <- sort(unique(site_labels))
    n_sites <- length(site_list)
    sigma_site_out <- matrix(NA, nrow = num_saved_iters, ncol = n_sites) # stores delta draws
  } else {
    site_list <- NULL
    n_sites <- 1
    sigma_site_out <- matrix(NA, nrow = num_saved_iters, ncol = 1)
  }
  
  # forests
  mu_forest  <- MakeForest(hypers_mu, opts_mu)
  tau_forest <- MakeForest(hypers_tau, opts_tau)
  
  # initials
  mu  <- rep(0, n)
  tau <- rep(0, n)
  sigma <- sd(Y)
  
  mu_forest$set_sigma(sigma)
  tau_forest$set_sigma(sigma)
  
  if (var_scaling) {
    delta_s <- rep(1.0, n_sites)
    v_s <- delta_s^2
    site_idx <- match(site_labels, site_list)
    site_weights <- 1 / v_s[site_idx]
  } else {
    delta_s <- 1.0
    v_s <- 1.0
    site_idx <- rep(1, n)
    site_weights <- rep(1, n)
  }
  
  alpha_delta <- alpha0
  beta_delta  <- beta0
  alpha_sigma <- alpha0
  beta_sigma  <- beta0
  
  save_index <- 1
  cat("Starting BARTharm sampling: num_iter =", num_iter, "thinning =", thinning_interval, "\n")
  
  for (iter in 1:num_iter) {
    ## mu update
    R_mu <- as.numeric(Y) - tau
    if (var_scaling) {
      mu <- mu_forest$do_gibbs_weighted(X_iqm_matrix, R_mu, site_weights, X_iqm_matrix, 1)
    } else {
      mu <- mu_forest$do_gibbs(X_iqm_matrix, R_mu, X_iqm_matrix, 1)
    }
    
    ## tau update
    R_tau <- as.numeric(Y) - mu
    if (var_scaling) {
      tau <- tau_forest$do_gibbs_weighted(X_bio_matrix, R_tau, site_weights, X_bio_matrix, 1)
    } else {
      tau <- tau_forest$do_gibbs(X_bio_matrix, R_tau, X_bio_matrix, 1)
    }
    
    ## residuals
    residuals <- as.numeric(Y) - mu - tau
    
    ## update v_s = delta^2 (no reference-site fix)
    if (var_scaling) {
      for (s_idx in seq_len(n_sites)) {
        mask <- (site_idx == s_idx)
        nj <- sum(mask)
        if (nj > 0) {
          res_j <- residuals[mask]
          shape_post <- alpha_delta + nj / 2
          scale_post <- beta_delta + sum(res_j^2) / (2 * sigma^2)
          v_draw <- 1 / rgamma(1, shape = shape_post, rate = scale_post)
          v_draw <- pmax(v_draw, 1e-12)
          v_s[s_idx] <- v_draw
          delta_s[s_idx] <- sqrt(v_draw)
        } else {
          # leave as-is
          v_s[s_idx] <- v_s[s_idx]
        }
      }
      site_weights <- 1 / v_s[ site_idx ]
    } else {
      v_s <- 1.0
      delta_s <- 1.0
      site_weights <- rep(1, n)
    }
    
    ## update sigma (global)
    if (var_scaling) {
      shape_sigma_post <- alpha_sigma + n / 2
      scale_sigma_post <- beta_sigma + 0.5 * sum( residuals^2 / v_s[ site_idx ] )
      sigma2_draw <- 1 / rgamma(1, shape = shape_sigma_post, rate = scale_sigma_post)
      sigma <- sqrt( pmax(sigma2_draw, 1e-12) )
    } else {
      shape_sigma_post <- alpha_sigma + n / 2
      scale_sigma_post <- beta_sigma + 0.5 * sum(residuals^2)
      sigma2_draw <- 1 / rgamma(1, shape = shape_sigma_post, rate = scale_sigma_post)
      sigma <- sqrt( pmax(sigma2_draw, 1e-12) )
    }
    
    ## push sigma into forests
    mu_forest$set_sigma(sigma)
    tau_forest$set_sigma(sigma)
    
    ## store
    if (iter %% thinning_interval == 0) {
      mu_out[save_index, ] <- mu
      tau_out[save_index, ] <- tau
      sigma_out[save_index] <- sigma
      if (var_scaling) {
        sigma_site_out[save_index, ] <- delta_s
      } else {
        sigma_site_out[save_index, 1] <- 1.0
      }
      save_index <- save_index + 1
    }
  } # end loop
  
  return(list(
    mu_out = mu_out,
    tau_out = tau_out,
    sigma_out = sigma_out,
    sigma_site_out = sigma_site_out,  # delta draws
    final = list(mu = mu, tau = tau, sigma = sigma, delta_s = delta_s, v_s = v_s)
  ))
}


# delta_draws: matrix draws x n_sites (sigma_site_out returned by function)
# sigma_draws: vector of length draws (sigma_out)
# returns list with aligned draws (same dims) and aligned sigma vector
align_deltas_posthoc <- function(delta_draws, sigma_draws) {
  n_draws <- nrow(delta_draws)
  n_sites <- ncol(delta_draws)
  delta_aligned <- matrix(NA, nrow = n_draws, ncol = n_sites)
  sigma_aligned <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    rowd <- delta_draws[d, ]
    gm <- exp(mean(log(rowd)))  # geometric mean for that draw
    delta_aligned[d, ] <- rowd / gm
    sigma_aligned[d] <- sigma_draws[d] * gm  # to preserve equivalence
  }
  list(delta_aligned = delta_aligned, sigma_aligned = sigma_aligned)
}


# Inputs:
# mu_out, tau_out: matrices draws x n_obs (returned by inference)
# delta_out: matrix draws x n_sites (either raw or aligned)
# sigma_out: vector draws length draws (either raw or aligned)
# Y_scaled: outcome used in fitting (scaled)
# mY, sY: original mean and sd used to scale Y (scalars)
# sites: integer vector mapping observation -> site index (1..n_sites)
# keep_idx: indices of draws to use (1:n_draws typically)

posterior_predict_bartharm <- function(mu_out, tau_out, delta_out, sigma_out,
                                       Y_scaled, mY, sY, sites, keep_idx = NULL) {
  if (is.null(keep_idx)) keep_idx <- seq_len(nrow(mu_out))
  n_draws <- length(keep_idx)
  n_obs <- ncol(mu_out)
  
  Yharm_draws_norm <- matrix(NA, nrow = n_draws, ncol = n_obs)
  
  for (i in seq_along(keep_idx)) {
    r <- keep_idx[i]
    mu_d  <- mu_out[r, ]
    tau_d <- tau_out[r, ]
    delta_site_d <- delta_out[r, ]    # per-site deltas
    delta_d <- delta_site_d[sites]
    
    # scaled-space harmonised draw
    Yharm_scaled_d <- (Y_scaled - mu_d - tau_d) / delta_d + tau_d
    
    # back-transform
    Yharm_draws_norm[i, ] <-  Yharm_scaled_d 
  }
  
  list(
    Yharm_draws_norm = Yharm_draws_norm,
    Yharm_pp_mean = colMeans(Yharm_draws_norm),
    Yharm_pp_median = apply(Yharm_draws_norm, 2, median)
  )
}


posterior_predict_bartharm_noscaling <- function(mu_out, tau_out, Y_scaled, sites, keep_idx = NULL) {
  if (is.null(keep_idx)) keep_idx <- seq_len(nrow(mu_out))
  n_draws <- length(keep_idx)
  n_obs <- ncol(mu_out)
  
  Yharm_draws_norm <- matrix(NA, nrow = n_draws, ncol = n_obs)
  
  for (i in seq_along(keep_idx)) {
    r <- keep_idx[i]
    mu_d  <- mu_out[r, ]
    tau_d <- tau_out[r, ]
    
    # scaled-space harmonised draw
    Yharm_scaled_d <- Y_scaled - mu_d 
    
    # back-transform
    Yharm_draws_norm[i, ] <-  Yharm_scaled_d 
  }
  
  list(
    Yharm_draws_norm = Yharm_draws_norm,
    Yharm_pp_mean = colMeans(Yharm_draws_norm),
    Yharm_pp_median = apply(Yharm_draws_norm, 2, median)
  )
}