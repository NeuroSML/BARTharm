# Loading libraries 
library(caret)
library(SoftBart)
library(dplyr)
library(matrixStats)
library(patchwork)
library(devtools)
library(tidyr)

# Loading functions
source("R/bartharm.R")
source("R/bartharm_inference.R")
source("R/get_data.R")
source("R/simulate_data.R")
source("R/normalise_data.R")
source("R/load_data.R")

# Define path to real dataset (must be an .RData file containing a data frame)
file_path <- 'data/real_data.RData'

# Directory where results and harmonized outputs will be saved
saving_path <- "examples/results/"

# Specify column names for biological covariates (e.g., age and sex)
bio_col <- c("Age", "Sex")

# Specify column names for image quality metrics (IQMs)
iqm_col <- c("snv", "cnr", "qi_1", "qi_2")

# Specify column names for the outcomes to be harmonized
outcomes_col <- c("NBV1", "NBV2")

# Specify the name of the subject ID column
id_col <- c("num_ID")

# Run BARTharm harmonization on real data
# - Loads and normalizes data
# - Runs BART-based Gibbs sampling to separate nuisance (mu) and signal (tau)
# - Saves harmonized outcomes and posterior draws

df_harmonised <- bartharm(
  simulate_data = FALSE,             # Use real data (not simulated)
  saving_path = saving_path,         # Output directory
  file_path = file_path,             # Input dataset path
  bio_col = bio_col,                 # Biological covariates
  iqm_col = iqm_col,                 # IQM covariates
  outcomes_col = outcomes_col,       # Outcome variables to harmonize
  id_col = id_col,                   # Subject ID column
  num_iter = 5000,                   # Number of MCMC iterations
  burn_in = 500,                     # Burn-in iterations
  thinning_interval = 2,             # Thinning
  num_tree_mu = 200,                 # Trees in mu forest (IQMs)
  num_tree_tau = 50,                 # Trees in tau forest (bio features)
  beta_mu = 2, beta_tau = 2,         # BART prior parameters
  gamma_mu = 0.95, gamma_tau = 0.95  # BART prior parameters
)

# Extract harmonized versions of the outcomes
harmonised_NBV1 <- df_harmonised$NBV1_harmonised
harmonised_NBV2 <- df_harmonised$NBV2_harmonised
