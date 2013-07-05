## prevalence

The prevalence package provides Frequentist and Bayesian methods useful in prevalence assessment studies.

To install or update, run:

    install.packages("prevalence")

IMPORTANT: the truePrev functions in the prevalence package call on JAGS (Just Another Gibbs Sampler), which therefore has to be available on the user's system. JAGS can be downloaded from http://mcmc-jags.sourceforge.net/.


## Development

To install the development version of the prevalence package, it is easiest to use the `devtools` package:

    install.packages("devtools")  # if needed..
    library(devtools)
    install_github("prevalence")