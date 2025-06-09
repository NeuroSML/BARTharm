# BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric

**BARTharm** is an R-based pipeline for harmonizing imaging-derived outcomes by separating biological signal from scanner-related effects using Image Quality Metrics (IQMs) instead of Scanner IDs. It uses Bayesian Additive Regression Trees (BART) with Gibbs sampling to estimate scanner components (`mu`, from IQMs) and biological effects (`tau`, from biological covariates).


### References

> Prevot E, et al., (2025). BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric. bioRxiv. Published online 2025. doi:10.1101/2025.06.04.657792 [link](https://www.biorxiv.org/content/10.1101/2025.06.04.657792v1)


---


## Installation

Install required R packages:

```r
install.packages(c("dplyr", "matrixStats", "caret"))
# If not already installed
devtools::install_github("theodds/SoftBart")  # for SoftBart
```

---

## Repository Structure

``` 
bartharm/
├── R/                        # Core R functions
│   ├── simulate_data.R       # Generate synthetic scanner-biased data
│   ├── load_data.R           # Load and process real data
│   ├── get_data.R            # Get data for BARTharm
│   ├── normalise_data.R      # Data normalization
│   ├── bartharm_inference.R  # Gibbs sampler for BARTharm
│   ├── bartharm.R            # Main harmonization function
├── data/                     # Real or example datasets
│   └── real_data.RData
├── examples/                 # Example usage scripts
│   ├── run_simulated.R
│   └── run_real.R
├── results/                  # Output harmonized datasets and posteriors
├── README.md                 # Project overview
└── .gitignore                # Ignore cache and intermediate files
```

---

## Usage

`bartharm()` is the main function used to perform harmonization via Bayesian Additive Regression Trees (BART). It supports both simulated and real datasets and outputs harmonized outcome variables by separating scanner-related nuisance variation (mu) from biological signal (tau).

The function:
- Loads and normalizes the data,
- Fits BART models to decompose nuisance vs. signal effects,
- Returns a data frame with both original and harmonized outcomes,
- Saves harmonized results and posterior samples to disk.

The `examples` directory contains code for running BARTharm harmonization on either simulated data or real data. To simply return the harmonized data, one can use the following:

```
df_harmonised <- bartharm(saving_path = saving_path, ... )
```

where `...` are the user-specified arguments needed for harmonization. The returned object df_harmonised contains:
- Original outcomes (e.g., NBV1, NBV2)
- Harmonized outcomes (e.g., NBV1_harmonised, NBV2_harmonised)

To extract them:

```
harmonised_NBV1 <- df_harmonised$NBV1_harmonised
harmonised_NBV2 <- df_harmonised$NBV2_harmonised
```



###  Automatic Saving to Disk

BARTharm automatically saves key outputs to the specified `saving_path` directory:

- Individual harmonized outcomes: For each outcome specified in outcomes_col, a separate .RData file is saved as `results/harmonised_<OutcomeName>.RData`
- Full harmonized dataset: The complete df_harmonised containing the original real/simulated data plus the harmonized outcomes is saved as `harmonised_simulated_df.RData` if using simulated data; `harmonised_realdata_df.RData` if using real data.
- Full Gibbs Chains: Posterior samples from the Gibbs sampler are saved as .RData files, including:
  - Mu chains (scanner-related nuisance effects) `results/mu_out_<OutcomeName>.RData`
  - Tau chains (biological signal effects) `results/tau_out_<OutcomeName>.RData`
  - Residual noise chains (posterior noise) `results/sigma_out_<OutcomeName>.RData`

These chains are crucial for diagnostics and uncertainty quantification. 

### Troubleshooting the outcome

Residual noise chains can be used to examine MCMC convergence and evaluate whether the chosen number of Gibbs samples (`num_iter`) or the burn in (`burn_in`) are adequate. You can plot the chain as follows 

```
load("results/sigma_out_<OutcomeName>.RData")
plot(sigma_out_<OutcomeName>, type = 'l', main = "Trace plot of residual noise for <OutcomeName>")
```



### How to Prepare the Real Data

To use the BARTharm harmonization pipeline, you must first prepare your dataset in a format compatible with the functions provided in this repository. Below are the key requirements and steps to ensure your data is correctly structured.

**Single Data Frame Format**
Ensure your dataset is a single .RData file containing one data frame. Rows correspond to observations and the colums correspond to your variables. This data frame must include all necessary variables:
- Biological covariates (e.g., Age, Sex)
- Image Quality Metrics (IQMs) (e.g., snv, cnr, qi_1, qi_2)
- Outcome variables to be harmonized (e.g., NBV1, NBV2)
- A unique subject identifier (e.g., num_ID)

**No Missing Data**
The dataset must be complete—all rows must have valid (non-missing) values for each variable of each observation.

**Consistent Column Naming**
The column names in your data frame must exactly match those you pass to the bartharm() function via the arguments.

**File Format**
Save your data frame as an .RData file using save() in R into the directory that you then pass to the `bartharm()` as the `file_path`. For example:
```
save(my_dataframe, file = "data/real_data.RData")
```

An example of a nicely formatted real dataset is provided in 'data/real_data.RData`.

### General recommendation

To ensure efficient and effective use of `bartharm()` across different datasets and applications, consider the following best practices for setting up your harmonization pipeline:

- Voxel-wise harmonization. If you intend to harmonize voxel-wise data (i.e., a large number of high-dimensional imaging features), it is strongly recommended to Sspecify one voxel/feature at a time in the `outcomes_col` argument and parallelize computation across voxels to reduce runtime.

- Harmonizing few summary features (e.g., IDPs). If you are working with a small number of imaging-derived phenotypes (IDPs) or summary metrics (e.g., NBV1, NBV2), you can safely specify all of them together in the outcomes_col list and will obtain results is feasible runtimes.

- Tuning BART Parameters. The following parameters control the flexibility and regularization of the BART priors used to model scanner-related nuisance effects (mu) and biological signal (tau):
  - num_tree_mu, num_tree_tau: The number of trees used in the BART ensemble for mu and tau, respectively. Increasing these increases model capacity and flexibility, but at the cost of higher computational burden. Use with caution in small datasets or low-signal settings.
  - beta_mu, beta_tau: Controls the variance of terminal node parameters; higher values shrink more toward zero. Increasing this reduces variance of th estimated effect, pushing towards homogenous effects. 
  - gamma_mu, gamma_tau: Controls the probability of splitting internal nodes; lower values lead to shallower trees (more regularization). Decreasing this, encourages shallower trees (i.e., more shrinkage), leading to more stable estimates in noisy or over-parameterized data.

    We recommend starting with the default values provided in this package, especially when he number of IQM covariates is larger than the number of biological covariates, and you have limited prior knowledge about appropriate levels of regularization.

