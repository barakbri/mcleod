# Monte Carlo methods for Estimating Latent OverDispersion of binomial and Poisson random variables.

The package contains several heirarchical Bayesian methods for estimating the mixing distribution of Binomial and Poisson sampels. Methods include regularized estimators of the CDF of the mixing distribution, confidence intervals for the CDF of the mixing distribution, and binomial and Poisson regression with a random intercept term following a general distribution.

# How to install
Install the package via:

```r
Sys.setenv(USE_CXX14 = 1)
# install.packages("devtools") # install only if devtools not previously installed
devtools::install_github('barakbri/mcleod', build_vignettes = TRUE)
```

# Where to start
After installation, type `browseVignettes(package = 'mcleod')` in the console and press the *PDF* link to see the package vigenette. The vignette goes over the statistical method, as well as simple and advanced use cases.


# Reproducing paper results
A github repository with scripts used to reproduce paper results is at: https://github.com/barakbri/Bayes_CI_Scripts

