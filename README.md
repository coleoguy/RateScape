# RateScape

RateScape detects and estimates branch-specific rate heterogeneity in discrete character evolution on phylogenies. It uses a spike-and-slab Bayesian mixture model with reversible-jump MCMC to identify lineages where the rate of discrete trait evolution differs from the background, and paints the results directly onto the tree.

## Installation

RateScape includes compiled C++ code (via Rcpp/RcppArmadillo), so you will need a working C++ compiler.

**On macOS:** Install Xcode Command Line Tools by running `xcode-select --install` in Terminal.

**On Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

**On Linux:** Install `r-base-dev` (Debian/Ubuntu) or the equivalent for your distribution.

Then install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("coleoguy/RateScape")
```

## Quick start

```r
library(RateScape)

# Simulate a tree and discrete character data with rate heterogeneity
tree <- ape::rcoal(50)
Q <- makeQ(k = 4)
sim <- simRateScape(tree, Q, lambda_sigma = c(0.8, 0.5),
                    pi = 0.8, return_rates = TRUE)

# Fit the spike-and-slab model
fit <- fitRateScape(tree, sim$data, Q = Q,
                    lambda_sigma = c(0.8, 0.5),
                    ngen = 5000, burn_in = 1000)

# Summarize and visualize
summarizeRates(fit)
plotRateTree(fit)
```

## Key functions

`makeQ()` builds rate matrices for symmetric Mk, ARD, or chromosome number models (gain, loss, polyploidy).

`simRateScape()` simulates discrete character evolution with optional branch-specific rate scalars drawn from a spike-and-slab prior.

`fitRateScape()` runs the MCMC sampler to estimate per-branch rate scalars, the mixing weight (pi), and slab variance (sigma-squared), returning Bayes factors and posterior samples.

`checkPrior()` visualizes the spike-and-slab prior to help with parameter specification.

`plotRateTree()` paints the phylogeny with estimated rate scalars for quick visual identification of rate-shifted lineages.

`rateCategories()` provides a maximum-likelihood alternative using discretized gamma rate categories.

## Citation

If you use RateScape in published work, please cite the accompanying paper:

Blackmon, H. RateScape: Bayesian detection of branch-specific rate variation in discrete character evolution. *Methods in Ecology and Evolution* (in prep).

## License

GPL (>= 2)
