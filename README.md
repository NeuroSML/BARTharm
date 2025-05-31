# BARTharm: MRI Harmonization Using Image Quality Metrics and Bayesian Non-parametric

**BARTharm** is an R-based pipeline for harmonizing imaging-derived outcomes by separating biological signal from scanner-related effects using Image Quality Metrics (IQMs) instead of Scanner IDs. It uses Bayesian Additive Regression Trees (BART) with Gibbs sampling to estimate scanner components (`mu`, from IQMs) and biological effects (`tau`, from biological covariates).

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


