
# This function normalizes the data and extracts the features to be harmonized.
# It applies quantile normalization to biological and IQM covariates while keeping the outcome variable.
#
# Arguments:
# - data_bio: Data frame of biological covariates
# - data_iqm: Data frame of IQM covariates
# - outcomes_col: Vector of column names corresponding to outcome variables
# - id_col: name of subject ID column

normalise_data <- function(data_bio, data_iqm, outcomes_col, id_col){
  
  # Apply quantile normalization to biological covariates, excluding num_ID and outcome columns
  norm_data_bio <- as.data.frame(quantile_normalize_bart(data_bio[, -which(names(data_bio) %in% c(id_col, outcomes_col))]))
  
  # Apply quantile normalization to IQM covariates, excluding num_ID
  norm_data_iqm <- as.data.frame(quantile_normalize_bart(data_iqm[, -which(names(data_iqm) %in% c(id_col))]))
  
  # Keep only the outcome variables and num_ID
  Y <- data_bio %>%
    dplyr::select(outcomes_col, id_col)
  
  # Return normalized data
  return(list("norm_data_bio" = norm_data_bio, "norm_data_iqm" = norm_data_iqm, "Y" = Y))
}