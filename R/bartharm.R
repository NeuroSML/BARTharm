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
# - site_col: Name of scanner/site ID column (required if var_scaling = TRUE).
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
# - var_scaling: Logical. If TRUE, harmonize variance across sites (requires site_col).

bartharm <- function(file_path = " ", saving_path = " ", save_format = "", simulate_data = TRUE, bio_col = c(), iqm_col = c(), outcomes_col = c(), id_col = c(), site_col = c(), n_subjects = 1000, linear_tau = TRUE, linear_mu = TRUE,
                     num_iter = 5000, burn_in = 500, thinning_interval = 2, num_tree_mu = 200, num_tree_tau = 50, beta_mu = 2, beta_tau = 2, gamma_mu = 0.95, gamma_tau = 0.95, var_scaling = FALSE){
  
  # Load or simulate data
  if(simulate_data){
    cat("Simulating data \n")
    data <- get_data(simulate = TRUE, saving_path = saving_path, save_format = save_format,  n_subjects = n_subjects, linear_tau = linear_tau, linear_mu = linear_mu)
    site_col <- data$site_col
    cat("Saved simulated data \n")
  }else{
    cat("Processing real data from", file_path ,"\n")
    data <- get_data(simulate = FALSE, filepath = file_path, save_format = save_format, saving_path = saving_path,  id_col = id_col, bio_col = bio_col, iqm_col = iqm_col, outcomes_col = outcomes_col, site_col = site_col, var_scaling = var_scaling)
  }
  
  # Extract normalized matrices and outcome
  X_bio_matrix <- data$X_bio_matrix
  X_iqm_matrix <- data$X_iqm_matrix
  Y <- as.data.frame(data$Y)
  Y_norm <- as.data.frame(data$Y_norm)
  original_means <- data$original_means
  original_sds <- data$original_sds
  df <- data$df
  
  if(var_scaling){
    cat("Using site information from column: ", site_col[1],"\n")
    sites <- X_iqm_matrix[, site_col[1]]
  } else{
    cat("Site information not available","\n")
    sites <- NULL
  }
  
  ll <- colnames(Y) # Names of outcome variables and ID column
  cat("Harmonizing: ", ll[1:(length(ll)-1)], "\n") # Skip ID column
  
  num_saved_iters <- floor(num_iter / thinning_interval)
  df_harmonised <- df  # Make a copy of the raw data to store harmonized results
  
  # Loop over each outcome variable (excluding ID column)
  for(i in 1:(length(ll)-1)){
    cat("Executing harmonization for feature: ", ll[i], "\n")
    
    # Set up BART hyperparameters for mu (scanner effect) and tau (biological effect)
    hypers_mu <- Hypers(X_iqm_matrix, Y_norm[,i], num_tree = num_tree_mu, beta = beta_mu, gamma = gamma_mu, normalize_Y = FALSE) 
    hypers_tau <- Hypers(X_bio_matrix, Y_norm[,i], num_tree = num_tree_tau, beta = beta_tau, gamma = gamma_tau, normalize_Y = FALSE)
    
    # Use default options for BART forests
    opts_mu <- Opts(update_sigma = FALSE)
    opts_tau <- Opts(update_sigma = FALSE)
    
    alpha0 = 0.01
    beta0 = 0.01
    
    # Run the Gibbs sampler for mu and tau
    bartharm_output <- bartharm_inference(num_iter, thinning_interval, X_iqm_matrix, X_bio_matrix, Y_norm[,i], hypers_mu, hypers_tau, opts_mu, opts_tau, var_scaling, sites, alpha0 = alpha0, beta0 = beta0)
    
    # Extract posterior draws
    mu_out <- bartharm_output$mu_out
    tau_out <- bartharm_output$tau_out
    sigma_out <- bartharm_output$sigma_out
    
    # Save posterior samples
    cat("Saving full posterior samples for feature: ", ll[i], "\n")
    saving_data(mu_out, file_name = paste0('mu_out_',ll[i]), saving_path, save_format = save_format)
    saving_data(tau_out, file_name = paste0('tau_out_',ll[i]), saving_path, save_format = save_format)
    saving_data(sigma_out, file_name = paste0('sigma_out_',ll[i]), saving_path, save_format = save_format)
    
    if(var_scaling){
      sigma_site_out <- bartharm_output$sigma_site_out
      save(file=paste0(saving_path, 'sigma_site_out_',ll[i],'.RData'), sigma_site_out)
    }
    
    # Compute posterior mean prediction
    y_pred <- colMeans(mu_out[(burn_in/thinning_interval+1):num_saved_iters, ]) + colMeans(tau_out[(burn_in/thinning_interval+1):num_saved_iters, ])
    
    # Evaluate RMSE between predicted and observed
    rmse_value <- rmse(Y_norm[,i], y_pred)
    cat("Prediction RMSE for feature ", ll[i], ": ", rmse_value, "\n")
    
    # Compute harmonized outcome by removing nuisance (mu) component
    cat("Evaluating harmonized feature: ", ll[i], "\n")
    
    if(var_scaling){
      
      print("Using variance scaling to calculated harmonized outcome \n")
      print("Aligning deltas \n")
      aligned <- align_deltas_posthoc(bartharm_output$sigma_site_out, bartharm_output$sigma_out)
      delta_aligned_draws <- aligned$delta_aligned
      sigma_aligned_draws <- aligned$sigma_aligned
      
      keep <- (burn_in/thinning_interval + 1):nrow(mu_out)
      print("Getting posterior mean and medians \n")
      pp <- posterior_predict_bartharm(mu_out = mu_out, tau_out = tau_out,
                                       delta_out = delta_aligned_draws, sigma_out = sigma_aligned_draws,
                                       Y_scaled = Y_norm[,i], mY = original_means[i], sY = original_sds[i],
                                       sites = sites, keep_idx = keep)
      
      
      Yharm_pp_mean_norm <- pp$Yharm_pp_mean
      Yharm_pp_median_norm <- pp$Yharm_pp_median

      
    } else{
      print("NOT Using variance scaling to calculated harmonized outcome \n")
      print("Getting posterior mean and medians \n")
      keep <- (burn_in/thinning_interval + 1):nrow(mu_out)
      pp <- posterior_predict_bartharm_noscaling(mu_out = mu_out, tau_out = tau_out,
                                       Y_scaled = Y_norm[,i], sites = sites, keep_idx = keep)
      
      Yharm_pp_mean_norm <- pp$Yharm_pp_mean
      Yharm_pp_median_norm <- pp$Yharm_pp_median
    }
    
    # Add harmonized and predicted values to the dataframe
    df_harmonised[, paste0(ll[i], "_harmonised_mean_raw")] <- Yharm_pp_mean_norm
    df_harmonised[, paste0(ll[i], "_harmonised_median_raw")] <- Yharm_pp_median_norm
    df_harmonised[, paste0(ll[i], "_predicted_raw")] <- y_pred
    
    y_harmonised_mean_original <- Yharm_pp_mean_norm * original_sds[i] + original_means[i]
    y_harmonised_median_original <- Yharm_pp_median_norm * original_sds[i] + original_means[i]
    y_pred_original <- y_pred * original_sds[i] + original_means[i]
    
    # Add harmonized and predicted values to the dataframe
    df_harmonised[, paste0(ll[i], "_harmonised_mean_original")] <- y_harmonised_mean_original
    df_harmonised[, paste0(ll[i], "_harmonised_median_original")] <- y_harmonised_median_original
    df_harmonised[, paste0(ll[i], "_predicted_original")] <- y_pred_original
    
    # Save harmonized outcome to disk
    cat("Saving harmonized feature at $harmonised_", ll[i] , " \n")
    saving_data(Yharm_pp_mean_norm, file_name = paste0('harmonised_mean_',ll[i], '_raw'), 
                saving_path, save_format = save_format)
    saving_data(Yharm_pp_median_norm, file_name = paste0('harmonised_median_',ll[i], '_raw'), 
                saving_path, save_format = save_format)
    saving_data(y_harmonised_mean_original, 
                file_name = paste0('harmonised_mean_',ll[i], '_original'), 
                saving_path, save_format = save_format)
    saving_data(y_harmonised_median_original, 
                file_name = paste0('harmonised_median_',ll[i], '_original'), 
                saving_path, save_format = save_format)
    
    cat("Saving predicted feature at $predicted_", ll[i] , " \n")
    saving_data(y_pred, file_name = paste0('predicted_',ll[i], '_raw'), 
                saving_path, save_format = save_format)
    saving_data(y_pred_original, 
                file_name = paste0('predicted_',ll[i], '_original'), 
                saving_path, save_format = save_format)
  }
  
  # Save the full harmonized dataframe (if running sequentially or simulating data)
  if(simulate_data){
    cat("Saving final harmonized dataset\n")
    saving_data(df_harmonised, "harmonised_simulated_df", saving_path, save_format = save_format)
  } else {
    if((length(ll)-1)>1){
      cat("Saving final harmonized dataset\n")
      saving_data(df_harmonised, "df_combined_harmonised_realdata", saving_path, save_format = save_format)
    }
  }
  
  return(df_harmonised)
  
}
