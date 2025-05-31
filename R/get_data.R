# This function retrieves and processes data for harmonization.
# If simulate = TRUE, it generates synthetic data using `simulate_data()`.
# If simulate = FALSE, it loads real data from a specified file using `load_data()`.
# In both cases, it applies quantile normalization to biological and IQM covariates using `normalise_data()`,
# and saves the intermediate datasets to disk.
#
# Arguments:
# - simulate: Logical. If TRUE, simulate synthetic data; otherwise, load real data.
# - filepath: Path to real data (.RData file) when simulate = FALSE.
# - saving_path: Directory to save intermediate and output files.
# - id_col: Name of subject ID column (character vector).
# - bio_col: Vector of biological covariate column names (only used if simulate = FALSE).
# - iqm_col: Vector of IQM covariate column names (only used if simulate = FALSE).
# - outcomes_col: Vector of outcome column names.
# - n_subjects: Number of subjects to simulate (default = 1000).
# - linear_tau: Logical. If TRUE, simulate outcome with linear biological effect.
# - linear_mu: Logical. If TRUE, model scanner effects linearly.


get_data <- function(simulate=FALSE, filepath = "", saving_path = "", id_col = c(), bio_col = c(), iqm_col = c(), outcomes_col = c(), 
                     n_subjects = 1000, linear_tau = TRUE, linear_mu = TRUE){
  if(simulate){
    data <- simulate_data(n_subjects, linear_tau, linear_mu)
    data_bio <- data$data_bio
    data_iqm <- data$data_iqm
    df <- data$simdata
    
    save(file=paste0(saving_path, 'simulated_df.RData'), df)
    
    id_col <- c("subid")
    outcomes_col <- c("outcome_simulated")
    
    cat("Normalising data \n")
    normalised_data <- normalise_data(data_bio, data_iqm, outcomes_col, id_col)
    
    norm_data_bio <- normalised_data$norm_data_bio
    norm_data_iqm <- normalised_data$norm_data_iqm
    
    save(file=paste0(saving_path, 'normalised_simdata_bio.RData'), norm_data_bio)
    save(file=paste0(saving_path, 'normalised_simdata_iqm.RData'), norm_data_iqm)
    
    
    
  }else{
    # Call function to load data
    data <- load_data(filepath, id_col, bio_col, iqm_col)
    data_bio <- data$data_bio
    data_iqm <- data$data_iqm
    df <- data$data
    
    save(file=paste0(saving_path, 'realdata_df.RData'), df)
    
    # Verify that the split happened correctly
    print("Bio column names:")
    print(colnames(data_bio))
    print("IQMs column names:")
    print(colnames(data_iqm))
    
    cat("Normalising data \n")
    normalised_data <- normalise_data(data_bio, data_iqm, outcomes_col, id_col)
    
    norm_data_bio <- normalised_data$norm_data_bio
    norm_data_iqm <- normalised_data$norm_data_iqm
    
    save(file=paste0(saving_path, 'normalised_realdata_bio.RData'), norm_data_bio)
    save(file=paste0(saving_path, 'normalised_realdata_iqm.RData'), norm_data_iqm)
    
  }
  
  return(list("X_bio_matrix" = as.matrix(norm_data_bio), "X_iqm_matrix" = as.matrix(norm_data_iqm), "Y" = normalised_data$Y, "df" = df))
}