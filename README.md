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

### Usage

`bartharm()` is the main function used to perform harmonization via Bayesian Additive Regression Trees (BART). It supports both simulated and real datasets and outputs harmonized outcome variables after removing scanner-related nuisance effects.

To simply return the harmonized data, one can use the following:

```
df_harmonised <- bartharm(saving_path = saving_path, ... )
```

where `...` are the user-specified arguments needed for harmonization. 

The `examples` directory contains code for running BARTharm harmonization on either simulated data or real data.  


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
