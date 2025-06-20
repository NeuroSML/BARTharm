# BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric

**BARTharm** is an R-based pipeline for harmonizing imaging-derived outcomes by separating biological signal from scanner-related effects using Image Quality Metrics (IQMs) instead of Scanner IDs. It uses Bayesian Additive Regression Trees (BART) with Gibbs sampling to estimate scanner components (`mu`, from IQMs) and biological effects (`tau`, from biological covariates).


### References

> Prevot E, et al., (2025). BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric. bioRxiv. Published online 2025. doi:10.1101/2025.06.04.657792 [link](https://www.biorxiv.org/content/10.1101/2025.06.04.657792v1)


---


## Installation

Install required R packages:

```r
# If not already installed
install.packages(c("dplyr", "matrixStats", "caret"))
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

## How to Prepare the Real Data

To use the BARTharm harmonization pipeline, you must first prepare your dataset in a format compatible with the functions provided in this repository. Below are the key requirements and steps to ensure your data is correctly structured.

**Extracting Image Quality Metrics (IQMs)**

Image Quality Metrics (IQMs) can be extracted using [MRIQC](https://mriqc.readthedocs.io/en/latest/), which is a tool to extract no-reference IQMs (image quality metrics) from structural (T1w and T2w), functional and diffusion MRI (magnetic resonance imaging) data. 

However, *BARTharm is not limited to a specific set of IQMs, or to a particular extraction tool such as MRIQC, or a particular imaging modality*. It can be flexibly applied to any collection of IQMs derived from any available quality assessment framework, provided they capture meaningful scanner-induced variability. Furthermore, BARTharm can accommodate any number or subset of metrics, the user just needs to specify the column names in the argument `iqm_col`.

**Single Data Frame Format**

Ensure your dataset is a single .RData file containing one data frame. Rows correspond to observations and the colums correspond to your variables. This data frame must include all necessary variables:
- Biological covariates (e.g., Age, Sex). 
- Image Quality Metrics (IQMs) (e.g., snv, cnr, qi_1, qi_2). 
- Outcome variables to be harmonized (e.g., NBV1, NBV2). 
- A unique subject identifier (e.g., num_ID). 

We recommend the User to use all available biological covariates and all available image quality covariates to perform BARTharm harmonization. 

**Consistent Column Naming**: The column names in your data frame must exactly match those you pass to the bartharm() function via the arguments `bio_col` for the biological covariates, `iqm_col` for the IQMs, `outcomes_col` for outcome variables to be harmonized, and `id_col` for the unique subject identifier (see `examples` directory).

**No Missing Data**: The dataset must be complete, all rows must have valid (non-missing) values for each variable of each observation.

**File Format**
You can provide your dataset in one of the following formats:

- .RData 
- .csv 
- .tsv 

Pass the path to one of the above file types as the `file_path` argument in `bartharm()`.

Examples of a correctly formatted real dataset are provided in the `data` directory.

## Usage

The `examples` directory contains code for running BARTharm harmonization on either simulated data or real data.

`bartharm()` is the main function used to perform harmonization via Bayesian Additive Regression Trees (BART). It supports both simulated and real datasets and outputs harmonized outcome variables by separating scanner-related nuisance variation (mu) from biological signal (tau). The simulated dataset which is generated in the `run_simulated.R` is the same described in the Simulation Framework 1 of the BARTharm paper. 

The `bartharm()` function:
- Loads/Simulates and normalizes the data,
- Fits BART models to decompose scanner vs. biological effects,
- Returns a data frame with both original and harmonized outcomes,
- Saves harmonized results and posterior samples to disk.

To simply return the harmonized data, one can use the following:

```
df_harmonised <- bartharm(saving_path = saving_path, ... )
```

where `...` are the user-specified arguments needed for harmonization and explained in `examples/run_simulated.R` or `examples/run_real.R`. 

The returned object df_harmonised is a data frame which contains:
- Original outcomes (e.g., NBV1, NBV2)
- Harmonized outcomes (e.g., NBV1_harmonised, NBV2_harmonised)
- Original real data or simulated data

To extract harmonized outcomes:

```
harmonised_NBV1 <- df_harmonised$NBV1_harmonised
harmonised_NBV2 <- df_harmonised$NBV2_harmonised
```

You can inspect the full returned dataframe with `head(df_harmonised)`.

###  Automatic Saving to Disk

The User can specify a format in the argument `save_format` between RData, csv and tsv for the saved dataframe. BARTharm automatically saves these key outputs to the specified `saving_path` directory:

- Original, pre-processed, Real Data or Simulated Data: The original dataset is saved as `simulated_df.<save_format>` for simulated data, `realdata_df.<save_format>` for real data, using the specified preferred `save_format`.
- Normalised biological covariate data and iqm covariate data. Original pre-processed, normalised, Real Data or Simulated Data divided into either bio or iqm are saved as `normalised_simdata_bio.<save_format>/ normalised_simdata_iqm.<save_format>` for simulated data, `normalised_realdata_bio.<save_format>/ normalised_realdata_iqm.<save_format>` for real data, using the specified preferred `save_format`.
- Individual harmonized outcomes: For each outcome specified in outcomes_col, a separate .RData file is saved as `results/harmonised_<OutcomeName>.RData`
- Full harmonized dataset: The complete df_harmonised containing the original data plus the harmonized outcomes is saved as `harmonised_simulated_df.<save_format>` for simulated data, `harmonised_realdata_df.<save_format>` for real data, using the specified preferred `save_format`.
- Full Gibbs Chains: Posterior samples from the Gibbs sampler are saved as .RData files, including:
  - Mu chains (scanner-related nuisance effects) `results/mu_out_<OutcomeName>.RData`
  - Tau chains (biological signal effects) `results/tau_out_<OutcomeName>.RData`
  - Residual noise chains (posterior noise) `results/sigma_out_<OutcomeName>.RData`

The `results` directory contains examplar outputs from running the `run_simulated.R` example script with `save_format = "RData`.


### Troubleshooting the outcome

Residual noise Gibbs chain can be used to examine MCMC convergence and evaluate whether the chosen number of Gibbs samples (`num_iter`) or the burn in (`burn_in`) are adequate. You can plot the chain as follows:


```
load("results/sigma_out_<OutcomeName>.RData")
plot(sigma_out_<OutcomeName>, type = 'l', main = "Trace plot of residual noise for <OutcomeName>")
```

Below is a simulated example of a well-mixing MCMC chain with a clear burn-in period. The chain starts off high and gradually stabilizes around its typical value after the first `~2000` iterations. The shaded region highlights the burn-in phase, which should be discarded before summarizing the results. This means that if you specified `burn_in > 2000` your estimates will be reliable. The specified `num_iter` is also high enough for the chain to reach convergence, resulting in a stable distribution of values.

![Alt text](examples/MCMC-burn-in-example.png)

If your chain instead looks like the one below — bouncing around slowly, stuck in places, or drifting a lot — this is a poorly mixing chain. It means the sampler is having trouble exploring the space properly. In this case, your results might not be trustworthy yet. When that happens, here’s what you can try

- Try increasing the number of iterations (`num_iter`) or letting it warm up longer by increasing the burn-in (`burn_in`).
- Try plotting again after the changes - you're looking for something that settles down and wiggles nicely like the first image.
- If that doesn't work, you can make it simpler by adjusting the settings to reduce the number of splits in the trees (lower `gamma`) or make the model more cautious (increase `beta`), or decrease number of trees.

![Alt text](examples/MCMC-poor-example.png)


### General recommendation

To ensure efficient and effective use of `bartharm()` across different datasets and applications, consider the following best practices for setting up your harmonization pipeline:

- **Voxel-wise harmonization**. If you intend to harmonize voxel-wise data (i.e., a large number of high-dimensional imaging features), it is strongly recommended to Sspecify one voxel/feature at a time in the `outcomes_col` argument and parallelize computation across voxels to reduce runtime.

- **Harmonizing few summary features (e.g., IDPs)**. If you are working with a small number of imaging-derived phenotypes (IDPs) or summary metrics (e.g., NBV1, NBV2), you can safely specify all of them together in the outcomes_col list and will obtain results is feasible runtimes.

- **Tuning BART Parameters**. The following parameters control the flexibility and regularization of the BART priors used to model scanner-related nuisance effects (mu) and biological signal (tau):
  - num_tree_mu, num_tree_tau: The number of trees used in the BART ensemble for mu and tau, respectively. Increasing these increases model capacity and flexibility, but at the cost of higher computational burden. Use with caution in small datasets or low-signal settings.
  - beta_mu, beta_tau: Controls the variance of terminal node parameters; higher values shrink more toward zero. Increasing this reduces variance of th estimated effect, pushing towards homogenous effects. 
  - gamma_mu, gamma_tau: Controls the probability of splitting internal nodes; lower values lead to shallower trees (more regularization). Decreasing this, encourages shallower trees (i.e., more shrinkage), leading to more stable estimates in noisy or over-parameterized data.

    We recommend starting with the default values provided in this package, especially when he number of IQM covariates is larger than the number of biological covariates, and you have limited prior knowledge about appropriate levels of regularization.


*For more guidance and explanations of the model we recommed the Users to look at the paper.*

