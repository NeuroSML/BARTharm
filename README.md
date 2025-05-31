# BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric

**BARTharm** is an R-based pipeline for harmonizing imaging-derived outcomes by separating biological signal from scanner-related effects using Image Quality Metrics (IQMs) instead of Scanner IDs. It uses Bayesian Additive Regression Trees (BART) with Gibbs sampling to estimate scanner components (`mu`, from IQMs) and biological effects (`tau`, from biological covariates).


### References

> Prevot, E., (2025). BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric.  [link]()


---


## Installation

Install required R packages:

```r
install.packages(c("dplyr", "matrixStats", "caret"))
# If not already installed
devtools::install_github("jaredsmurray/SoftBart")  # for SoftBart
```

---

## Repository Structure

``` 
bartharm/
├── R/                      # Core R functions
│   ├── simulate_data.R     # Generate synthetic scanner-biased data
│   ├── load_data.R         # Load and split real data
│   ├── normalise_data.R    # Quantile normalization
│   ├── bartharm_inference.R # Gibbs sampler for BARTharm
│   ├── bartharm.R          # High-level harmonization function
│   └── utils.R             # Utilities (e.g., rmse)
├── data/                   # Real or example datasets
│   └── real_data.RData
├── examples/               # Example usage scripts
│   ├── run_simulated.R
│   └── run_real.R
├── results/                # Output harmonized datasets and posteriors
├── README.md               # Project overview
└── .gitignore              # Ignore cache and intermediate files
```

---

### Usage

`bartharm()` is the main function used to perform harmonization via Bayesian Additive Regression Trees (BART). It supports both simulated and real datasets and outputs harmonized outcome variables after removing scanner-related nuisance effects.

To simply return the harmonized data, one can use the following:

```
df_harmonised <- bartharm(saving_path = saving_path, ... )
```

where `...` are the user-specified arguments needed for harmonization. 

The `examples` directory contains code for running BARTharm harmonization on either simulated data or real data.  






