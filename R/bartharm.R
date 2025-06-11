# This function performs BART-based harmonization on either simulated or real data.
# It prepares the data (via `get_data()`), then applies the BARTharm algorithm using Gibbs sampling
# to separate and remove scanner-related nuisance variation (mu) from biological signal (tau).
# Posterior samples and harmonized outcomes are saved to disk.
# Arguments:
# - file_path: Path to real data file (.RData), used when simulate_data = FALSE.
# - saving_path: Directory path to save intermediate and final outputs.
# - simulate_data: Logical. If TRUE, generates and harmonizes simulated data.
# - bio_col: Vector of biological covariate column names (for real data).
# - iqm_col: Vector of IQM covariate column names (for real data).
# - outcomes_col: Vector of outcome column names.
# - id_col: Name of subject ID column.
# - n_subjects: Number of simulated subjects (default = 1000).
# - linear_tau: Logical. If TRUE, simulate outcome with linear biological effects.
# - linear_mu: Logical. If TRUE, simulate scanner effects linearly.
# - num_iter: Total number of MCMC iterations.
# - burn_in: Number of iterations to discard as burn-in.
# - thinning_interval: Interval for saving MCMC samples.
# - num_tree_mu: Number of trees in the mu forest (IQMs).
# - num_tree_tau: Number of trees in the tau forest (biological covariates).
# - beta_mu, beta_tau: BART prior parameters for mu and tau forests.
# - gamma_mu, gamma_tau: BART prior parameters controlling sparsity.


bartharm <- function(file_path = " ", saving_path = " ", save_format = "", simulate_data = TRUE, bio_col = c(), iqm_col = c(), outcomes_col = c(), id_col = c(), n_subjects = 1000, linear_tau = TRUE, linear_mu = TRUE,
                     num_iter = 5000, burn_in = 500, thinning_interval = 2, num_tree_mu = 200, num_tree_tau = 50, beta_mu = 2, beta_tau = 2, gamma_mu = 0.95, gamma_tau = 0.95){
  
  # Load or simulate data
  if(simulate_data){
    cat("Simulating data \n")
    data <- get_data(simulate = TRUE, saving_path = saving_path, save_format = save_format,  n_subjects = n_subjects, linear_tau = linear_tau, linear_mu = linear_mu)
    cat("Saved simulated data \n")
  }else{
    cat("Processing real data from", file_path ,"\n")
    data <- get_data(simulate = FALSE, filepath = file_path, save_format = save_format, saving_path = saving_path,  id_col = id_col, bio_col = bio_col, iqm_col = iqm_col, outcomes_col = outcomes_col)
  }
  
  # Extract normalized matrices and outcome
  X_bio_matrix <- data$X_bio_matrix
  X_iqm_matrix <- data$X_iqm_matrix
  Y <- as.data.frame(data$Y)
  df <- data$df
  
  ll <- colnames(Y) # Names of outcome variables and ID column
  cat("Harmonizing: ", ll[1:(length(ll)-1)], "\n") # Skip ID column
  
  num_saved_iters <- ceiling(num_iter / thinning_interval) # Number of posterior samples saved
  
  df_harmonised <- df  # Make a copy of the raw data to store harmonized results
  
  # Loop over each outcome variable (excluding ID column)
  for(i in 1:(length(ll)-1)){
    cat("Executing harmonization for feature: ", ll[i], "\n")
    
    # Set up BART hyperparameters for mu (scanner effect) and tau (biological effect)
    hypers_mu <- Hypers(X_iqm_matrix, Y[,i], num_tree = num_tree_mu, beta = beta_mu, gamma = gamma_mu) 
    hypers_tau <- Hypers(X_bio_matrix, Y[,i], num_tree = num_tree_tau, beta = beta_tau, gamma = gamma_tau)
    
    # Use default options for BART forests
    opts_mu <- Opts()
    opts_tau <- Opts()
    
    # Run the Gibbs sampler for mu and tau
    bartharm_output <- bartharm_inference(num_iter, thinning_interval, X_iqm_matrix, X_bio_matrix, Y[,i], hypers_mu, hypers_tau, opts_mu, opts_tau)
    
    # Extract posterior draws
    mu_out <- bartharm_output$mu_out
    tau_out <- bartharm_output$tau_out
    sigma_out <- bartharm_output$sigma_out
    
    # Save posterior samples
    cat("Saving full posterior samples for feature: ", ll[i], "\n")
    
    save(file=paste0(saving_path, 'mu_out_',ll[i],'.RData'), mu_out)
    save(file=paste0(saving_path, 'tau_out_',ll[i],'.RData'), tau_out)
    save(file=paste0(saving_path, 'sigma_out_',ll[i],'.RData'), sigma_out)
    
    # Compute posterior mean prediction
    y_pred <- colMeans(mu_out[(burn_in:num_saved_iters), ]) + colMeans(tau_out[(burn_in:num_saved_iters), ])
    
    # Evaluate RMSE between predicted and observed
    rmse_value <- rmse(Y[,i], y_pred)
    cat("Prediction RMSE for feature ", ll[i], ": ", rmse_value, "\n")
    
    # Compute harmonized outcome by removing nuisance (mu) component
    cat("Evaluating harmonized feature: ", ll[i], "\n")
    y_harmonised <- Y[,i] - colMeans(mu_out[(burn_in:num_saved_iters), ])
    
    # Add harmonized and predicted values to the dataframe
    df_harmonised[, paste0(ll[i], "_harmonised")] <- y_harmonised
    df_harmonised[, paste0(ll[i], "_predicted")] <- y_pred
    
    # Save harmonized outcome to disk
    cat("Saving harmonized feature at $harmonised_", ll[i] , " \n")
    save(file=paste0(saving_path, 'harmonised_',ll[i],'.RData'), y_harmonised)
  }
  
  # Save the full harmonized dataframe
  cat("Saving final harmonized dataset\n")
  if(simulate_data){
    #save(file=paste0(saving_path, 'harmonised_simulated_df.RData'), df_harmonised)
    saving_data(df_harmonised, "harmonised_simulated_df", saving_path, save_format = save_format)
  }else{
    #save(file=paste0(saving_path, 'harmonised_realdata_df.RData'), df_harmonised)
    saving_data(df_harmonised, "harmonised_realdata_df", saving_path, save_format = save_format)
  }
  
  return(df_harmonised)
  
}
