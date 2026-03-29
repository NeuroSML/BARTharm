# BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric

**BARTharm** is an R-based pipeline for harmonizing imaging-derived outcomes by separating biological signal from scanner-related effects using Image Quality Metrics (IQMs) instead of Scanner IDs. It uses Bayesian Additive Regression Trees (BART) with Gibbs sampling to estimate scanner components (`mu`, from IQMs) and biological effects (`tau`, from biological covariates).


### References

> Prevot E, et al., (2025). BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric. bioRxiv. Published online 2025. doi:10.1101/2025.06.04.657792 [link](https://www.biorxiv.org/content/10.1101/2025.06.04.657792v1)


---


## Installation

Install required R packages:

```r
# If not already installed
install.packages(c("dplyr", "matrixStats", "caret", "MCMCpack"))
devtools::install_github("theodds/SoftBart")  # for SoftBart
```

---

## Repository Structure

``` 
bartharm/
├── R/                                   # Core R functions
│   ├── simulate_data.R                  # Generate synthetic scanner-biased data
│   ├── load_data.R                      # Load and process real data
│   ├── get_data.R                       # Get data for BARTharm
│   ├── normalise_data.R                 # Data normalization
│   ├── bartharm_inference.R             # Gibbs sampler for BARTharm
│   ├── bartharm.R                       # Main harmonization function
|   |__ combine_harmonised_outcomes.R    # Postprocessing function 
├── data/                                # Real or example datasets
│   └── real_data.RData
│   └── real_data.csv
│   └── real_data.tsv
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
- If available, a unique scanner ID.

We recommend the User to use all available biological covariates and all available image quality covariates to perform BARTharm harmonization. 

**Consistent Column Naming**: The column names in your data frame must exactly match those you pass to the bartharm() function via the arguments `bio_col` for the biological covariates, `iqm_col` for the IQMs, `outcomes_col` for outcome variables to be harmonized, `id_col` for the unique subject identifier, and `site_col` for the scanner/site ID (if available) (see `examples` directory).

**No Missing Data**: The dataset must be complete, all rows must have valid (non-missing) values for each variable of each observation.

**File Format**
You can provide your dataset in one of the following formats:

- .RData 
- .csv 
- .tsv 
- .rds

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

If `var_scaling` is set to TRUE, user needs to also provide scanner IDs and the model harmonises both the mean and the variance. Otherwise just the mean.

The returned object df_harmonised is a data frame which contains:
- Original outcomes 
- Harmonized outcomes 
- Predicted outcomes 
- Original real data or simulated data

You can inspect the full returned dataframe with `head(df_harmonised)`.

On real data, when running `bartharm()` sequentially, i.e., `length(outcomes_col) > 1`, the returned `df_harmonised` will already include all harmonised outcomes provided in `outcomes_col` and is saved to disk at `<saving_path>/df_combined_harmonised_realdata.<save_format>`.
When instead running `bartharm()` in parallel, one outcome at a time, `df_harmonised` will only include the current outcome being harmonised. 

###  Automatic Saving to Disk

The User can specify a format in the argument `save_format` between RData, csv, rds and tsv for the saved dataframe. BARTharm automatically saves these key outputs to the specified `saving_path` directory:

- Original, pre-processed, Real Data or Simulated Data: The original dataset is saved as `<saving_path>/simulated_df.<save_format>` for simulated data, `<saving_path>/filtered_realdata_df.<save_format>` for real data, using the specified preferred `save_format`.
- Normalised biological covariate data and iqm covariate data. Original pre-processed, normalised, Real Data or Simulated Data divided into either bio or iqm are saved as `<saving_path>/normalised_simdata_bio.<save_format>/ normalised_simdata_iqm.<save_format>` for simulated data, `<saving_path>/normalised_realdata_bio.<save_format>/ <saving_path>/normalised_realdata_iqm.<save_format>` for real data, using the specified preferred `save_format`.
- Individual harmonized and predicted outcomes: For each outcome specified in outcomes_col, separate files are saved as `<saving_path>/harmonised_<OutcomeName>_raw.<save_format>` or `<saving_path>/harmonised_<OutcomeName>_original.<save_format>`. `_raw` represents the harmonised outcome in the 0-1 normalised scale, and should be used if doing meta-analysis between different harmonised datasets. `_original` represents the harmonised outcome in the original outcome scale and can be compared to the observed outcome to assess harmonization. If `var_scaling = TRUE`, `<saving_path>/harmonised_<OutcomeName>_scaled_raw.<save_format>` or `<saving_path>/harmonised_<OutcomeName>_scaled_original.<save_format>` are also saved, where variance scaling is applied. The predicted outcomes are also saved in both scaled, which can be used to assess the goodness of the BARTharm fit to the data.
- Full Gibbs Chains: Posterior samples from the Gibbs sampler are saved in the pre-specified format, including:
  - Mu chains (scanner-related nuisance effects) `<saving_path>/mu_out_<OutcomeName>.<save_format>`
  - Tau chains (biological signal effects) `<saving_path>/tau_out_<OutcomeName>.<save_format>`
  - Residual noise chains (posterior noise) `<saving_path>/sigma_out_<OutcomeName>.<save_format>`
  - If `var_scaling = TRUE`, scanner/site-specific variances are also saved `<saving_path>/sigma_site_out_<OutcomeName>.<save_format>`
- Full harmonized dataset: If `bartharm()` is run sequentially, the complete df_harmonised containing the original data plus the harmonized and predicted outcomes is saved as `<saving_path>/df_combined_harmonised_realdata.<save_format>` for real data, `<saving_path>/harmonised_simulated_df.<save_format>` for simulated data, using the specified preferred `save_format`. Otherwise, if `bartharm()` is run in parallel, one outcome at a time, we provide a post-processing script to combined everything into a single dataframe once all outcomes have been harmonised. To use it, make sure all output files are correctly saved into the pre-specified `saving_path` and then run 
```
df_harmonised <- combine_harmonised_outcomes(saving_path, save_format)
```



### Troubleshooting the outcome

Residual noise Gibbs chain can be used to examine MCMC convergence and evaluate whether the chosen number of Gibbs samples (`num_iter`) or the burn in (`burn_in`) are adequate. You can plot the chain as follows:


```
load("<saving_path>/sigma_out_<OutcomeName>.<save_format>")
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

- **Voxel-wise harmonization**. If you intend to harmonize voxel-wise data (i.e., a large number of high-dimensional imaging features), it is strongly recommended to specify one voxel/feature at a time in the `outcomes_col` argument and parallelize computation across voxels to reduce runtime.

- **Harmonizing few summary features (e.g., IDPs)**. If you are working with a small number of imaging-derived phenotypes (IDPs) or summary metrics (e.g., NBV1, NBV2), you can safely specify all of them together in the outcomes_col list and will obtain results is feasible runtimes.

- **Tuning BART Parameters**. The following parameters control the flexibility and regularization of the BART priors used to model scanner-related nuisance effects (mu) and biological signal (tau):
  - num_tree_mu, num_tree_tau: The number of trees used in the BART ensemble for mu and tau, respectively. Increasing these increases model capacity and flexibility, but at the cost of higher computational burden. Use with caution in small datasets or low-signal settings.
  - beta_mu, beta_tau: Controls the variance of terminal node parameters; higher values shrink more toward zero. Increasing this reduces variance of th estimated effect, pushing towards homogenous effects. 
  - gamma_mu, gamma_tau: Controls the probability of splitting internal nodes; lower values lead to shallower trees (more regularization). Decreasing this, encourages shallower trees (i.e., more shrinkage), leading to more stable estimates in noisy or over-parameterized data.

    We recommend starting with the default values provided in this package, especially when he number of IQM covariates is larger than the number of biological covariates, and you have limited prior knowledge about appropriate levels of regularization.


*For more guidance and explanations of the model we recommed the Users to look at the paper.*

