
# This function loads the data and splits it into biological covariates and IQMs covariates
# Data processing should be done beforehand to ensure all columns contain numerical variables or factors.
#
# Arguments:
# - filepath: Path to the dataset
# - id_col: Column name containing subject IDs
# - bio_col: Vector of column names corresponding to biological covariates
# - iqm_col: Vector of column names corresponding to IQM covariates

load_data <- function(filepath, id_col, bio_col, iqm_col){
  
  # Load the dataset and remove rows with missing values
  loaded_data <- load(file=filepath)
  df <- get(loaded_data)
  df <- na.omit(df)
  
  # Select biological covariate columns along with numerical ID
  data_bio <-  df %>%
    dplyr::select(id_col, bio_col, outcomes_col)
  
  # Select IQM covariate columns along with numerical ID
  data_iqm <-  df %>%
    dplyr::select(id_col, iqm_col)
  
  # Return the processed datasets
  return(list("data" = df, "data_bio" = data_bio, "data_iqm" = data_iqm))
}