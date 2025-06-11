
# This function loads the data and splits it into biological covariates and IQMs covariates
# Data processing should be done beforehand to ensure all columns contain numerical variables or factors.
#
# Arguments:
# - filepath: Path to the dataset
# - id_col: Column name containing subject IDs
# - bio_col: Vector of column names corresponding to biological covariates
# - iqm_col: Vector of column names corresponding to IQM covariates

load_data <- function(filepath, id_col, bio_col, iqm_col){

  # Determine file extension
  ext <- tools::file_ext(filepath)
  
  # Load the dataset based on file type
  if (ext == "RData" || ext == "rda") {
    cat("Loading RData file\n")
    loaded_name <- load(file = filepath)
    df <- get(loaded_name)
  } else if (ext == "csv") {
    cat("Loading csv file\n")
    df <- read.csv(filepath, stringsAsFactors = TRUE)
  } else if (ext == "tsv") {
    cat("Loading tsv file\n")
    df <- read.delim(filepath, stringsAsFactors = TRUE)
  } else {
    stop("Unsupported file type: must be .RData, .csv, or .tsv")
  }

  # Remove rows with missing values
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