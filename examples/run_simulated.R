# Loading libraries 
library(caret)
library(SoftBart)
library(dplyr)
library(matrixStats)

# Loading functions
source("R/bartharm.R")
source("R/bartharm_inference.R")
source("R/get_data.R")
source("R/simulate_data.R")
source("R/normalise_data.R")
source("R/load_data.R")
source("R/saving_data.R")


# Define the directory where results and intermediate files will be saved
saving_path <- "results/"

# Specify the saving format
save_format = "RData" # can also be csv or tsv

# Run the BARTharm pipeline on simulated data
# - Simulates 1000 subjects with linear outcome and scanner effects
# - Uses 5000 MCMC iterations with a burn-in of 500 and thinning interval of 2
# - Sets number of trees and BART priors for both mu (scanner-related) and tau (biological) forests
df_harmonised <- bartharm(
  simulate_data = TRUE,              # Use simulated data
  saving_path = saving_path,         # Save results to this path
  save_format = save_format,         # Saving format
  n_subjects = 1000,                 # Number of subjects to simulate
  linear_tau = TRUE,                 # Linear/Non-Linear biological effects
  linear_mu = TRUE,                  # Linear/Non-Linear scanner effects
  num_iter = 500,                   # Total MCMC iterations
  burn_in = 50,                     # Number of burn-in samples to discard
  thinning_interval = 2,             # Thinning
  num_tree_mu = 200,                 # Trees in mu forest (IQMs)
  num_tree_tau = 200,                # Trees in tau forest (biological)
  beta_mu = 2, beta_tau = 2,         # BART prior parameters
  gamma_mu = 0.95, gamma_tau = 0.95  # BART prior parameters
)

# Extract the harmonized outcome variable from the output
harmonised_outcome <- df_harmonised$outcome_simulated_harmonised

# To evaluate the goodness of harmonization you can compare harmonised_outcome with the simulated clean data df_harmonised$outcome_clean