### The 'prevalence' package

The prevalence package provides Frequentist and Bayesian methods useful in prevalence assessment studies. Several methods are available for estimating True Prevalence (TP) from Apparent Prevalence (AP).

Available functions:

* `propCI`: Derive confidence intervals for a prevalence estimate
* `truePrev`: Estimate TP from AP obtained by testing individual samples with a single test
* `truePrevMulti`: Estimate TP from AP obtained by testing individual samples with multiple tests, using a conditional probability scheme
* `truePrevMulti2`: Estimate TP from AP obtained by testing individual samples with multiple tests, using a covariance scheme
* `truePrevPools`: Estimate TP from AP obtained by testing pooled samples
* `betaPERT`: Calculate the parameters of a Beta-PERT distribution
* `betaExpert`: Calculate the parameters of a Beta distribution based on expert opinion 

To install..

* .. the latest released version: `install.packages("prevalence")`
* .. the latest development version: `install_github("brechtdv/prevalence")`

IMPORTANT: the truePrev functions in the prevalence package call on JAGS (Just Another Gibbs Sampler), which therefore has to be available on the user's system. JAGS can be downloaded from http://mcmc-jags.sourceforge.net/.

Function `truePrev` is also available as an online Shiny application: http://users.ugent.be/~bdvleess/R/prevalence/shiny/

More information and tutorials are available at http://users.ugent.be/~bdvleess/R/prevalence/