
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

bartharm_inference <- function(num_iter, thinning_interval, X_iqm_matrix, X_bio_matrix, Y, hypers_mu, hypers_tau, opts_mu, opts_tau){
  
  num_saved_iters <- ceiling(num_iter / thinning_interval) # Number of posterior samples saved
  
  # Initialize matrices to store posterior samples
  mu_out <- matrix(NA, nrow = num_saved_iters, ncol = nrow(X_iqm_matrix))
  tau_out  <- matrix(NA, nrow = num_saved_iters, ncol = nrow(X_bio_matrix))
  sigma_out <- numeric(num_saved_iters)
  
  # Create BART forest objects for mu and tau estimation
  mu_forest <- MakeForest(hypers_mu, opts_mu)
  tau_forest <- MakeForest(hypers_tau, opts_tau)
  
  # Initialize tau using a linear model fit (alternatively this can just be initialised as zero)
  tau_lm <- lm(as.numeric(Y) ~ X_bio_matrix)
  tau <- as.numeric(predict(tau_lm, newdata = as.data.frame(X_bio_matrix)))
  
  save_index <- 1  # Index for storing samples
  
  cat("Starting sampling\n")
  
  # Main Gibbs sampling loop
  for(iter in 1:num_iter) {
    
    R <- (as.numeric(Y) - tau)  # Update residuals for mu estimation
    
    mu <- mu_forest$do_gibbs(X_iqm_matrix, R, X_iqm_matrix, 1)  # Update mu using Gibbs sampling
    
    sigma <- mu_forest$get_sigma()  # Get current estimate of sigma
    tau_forest$set_sigma(sigma)  # Set sigma in tau forest
    
    R <- (as.numeric(Y) - mu)  # Update residuals for tau estimation
    
    tau <- tau_forest$do_gibbs(X_bio_matrix, R, X_bio_matrix, 1)  # Update tau using Gibbs sampling
    
    # Save samples at thinning intervals
    if (iter %% thinning_interval == 0) {
      mu_out[save_index,] <- mu
      tau_out[save_index,]  <- tau
      sigma_out[save_index]  <- sigma
      save_index <- save_index + 1 
    }
  }
  
  # Return estimated parameters
  return(list("mu_out" = mu_out, "tau_out" = tau_out, "sigma_out" = sigma_out))
}