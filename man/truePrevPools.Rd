\name{truePrevPools}
\alias{truePrevPools}

\title{Estimate true prevalence from pooled samples}

\description{
Bayesian estimation of true prevalence from apparent prevalence obtained by
testing \emph{pooled} samples.
}

\usage{
truePrevPools(x, n, SE = 1, SP = 1, prior = c(1, 1),
              nchains = 2, burnin = 10000, update = 10000,
              verbose = FALSE)
}

\arguments{
  \item{x}{The vector of indicator variables, indicating whether a pool was
           positive (\code{"1"}) or negative (\code{"0"})}
  \item{n}{The vector of pool sizes}
  \item{SE, SP}{The prior distribution for sensitivity (SE) and specificity (SP); see 'Details' below for specification of these distributions}
  \item{prior}{The parameters of the prior Beta distribution for true prevalence; defaults to \code{c(1, 1)}}
  \item{nchains}{The number of chains used in the estimation process; \code{nchains} must be \eqn{\ge 2}}
  \item{burnin}{The number of discarded model iterations; defaults to 10,000}
  \item{update}{The number of withheld model iterations; defaults to 10,000}
  \item{verbose}{Logical flag, indicating if JAGS process output should be printed to the R console; defaults to \code{FALSE}}
}

\details{
  \code{truePrevPools} calls on \pkg{JAGS}/\pkg{\link[rjags:rjags]{rjags}} 
  to estimate the true prevalence from the apparent prevalence in a Bayesian
  framework. The default model, in BUGS language, is given below. To see the
  actual fitted model, see the model slot of the
  \link[prevalence:prev-class]{prev}-object.\cr
  \preformatted{
  model {
    for (i in 1:N) {
      x[i] ~ dbern(AP[i])
      AP[i] <- SEpool[i] * (1 - pow(1 - TP, n[i])) + (1 - SPpool[i]) * pow(1 - TP, n[i])
      SEpool[i] <- 1 - (pow(1 - SE, n[i] * TP) * pow(SP, n[i] * (1 - TP)))
      SPpool[i] <- pow(SP, n[i])
    }
  # SE ~ user-defined (see below)
  # SP ~ user-defined (see below)
  TP ~ dbeta(prior[1], prior[2])
  }
  }
  
  The test sensitivity (\code{SE}) and specificity (\code{SP}) can be specified
  by the user, independently, as one of \code{"fixed"}, \code{"uniform"},
  \code{"beta"}, \code{"pert"}, or \code{"beta-expert"}, with \code{"fixed"} as
  the default. Note that \code{SE} and \code{SP} must correspond to the test
  characteristics for testing individual samples; \code{truePrevPools} will
  calculate \code{SEpool} and \code{SPpool}, the sensitivity and specificitiy
  for testing pooled samples, based on Boelaert et al. (2000).

  Distribution parameters can be specified in a \emph{named} \code{list()}
  as follows:

\itemize{
  \item{\strong{Fixed:  }}{\code{list(dist = "fixed", par)}}
  \item{\strong{Uniform:  }}{\code{list(dist = "uniform", min, max)}}
  \item{\strong{Beta:  }}{\code{list(dist = "beta", alpha, beta)}}
  
  \item{\strong{PERT:  }}{\code{list(dist = "pert", method, a, m, b, k)}\cr
    \code{'method'} must be \code{"classic"} or \code{"vose"};\cr
    \code{'a'} denotes the pessimistic (minimum) estimate, \code{'m'} the most likely estimate, and \code{'b'} the optimistic (maximum) estimate;\cr
    \code{'k'} denotes the scale parameter.\cr
    See \code{\link{betaPERT}} for more information on Beta-PERT parametrization.}
    
  \item{\strong{Beta-Expert:  }}{\code{list(dist = "beta-expert", mode, mean, lower, upper, p)}\cr
    \code{'mode'} denotes the most likely estimate, \code{'mean'} the mean estimate;\cr
    \code{'lower'} denotes the lower bound, \code{'upper'} the upper bound;\cr
    \code{'p'} denotes the confidence level of the expert.\cr
	Only \code{mode} or \code{mean} should be specified; \code{lower} and
	\code{upper} can be specified together or alone.\cr
    See \code{\link{betaExpert}} for more information on Beta-Expert parameterization.}
}

  For Uniform, Beta and Beta-PERT distributions, BUGS-style short-hand notation
  is also allowed:
  
\itemize{
  \item{\strong{Uniform: }}{\code{~dunif(min, max)}}
  \item{\strong{Beta: }}{\code{~dbeta(alpha, beta)}}
  \item{\strong{Beta-PERT: }}{\code{~dpert(min, mode, max)}}
}}

\value{
An object of class \link[prevalence:prev-class]{prev}.
}

\note{
Markov chain Monte Carlo sampling in \code{truePrevPools} is performed by
\pkg{JAGS (Just Another Gibbs Sampler)} through the
\pkg{\link[rjags:rjags]{rjags}} package. JAGS can be downloaded from
\url{https://mcmc-jags.sourceforge.io/}.
}

\author{
Brecht Devleesschauwer <\email{brechtdv@gmail.com}>
}

\references{
\itemize{
  \item{Speybroeck N, Williams CJ, Lafia KB, Devleesschauwer B, Berkvens D (2012) Estimating the prevalence of infections in vector populations using pools of samples. \emph{Med Vet Entomol} \strong{26}:361-371}
  \item{Boelaert F, Walravens K, Biront P, Vermeersch JP, Berkvens D, Godfroid J (2000) Prevalence of paratuberculosis (Johne's disease) in the Belgian cattle population. \emph{Vet Microbiol} \strong{77}:269-281}
}}

\seealso{
  \pkg{\link[coda:mcmc]{coda}} for various functions that can be applied to the \code{prev@mcmc} object\cr
  \code{\link{truePrev}}: estimate true prevalence from apparent prevalence obtained by testing \emph{individual} samples with a single test\cr
  \code{\link{truePrevMulti}}: estimate true prevalence from apparent prevalence obtained by testing \emph{individual} samples with multiple tests, using a conditional probability scheme\cr
  \code{\link{truePrevMulti2}}: estimate true prevalence from apparent prevalence obtained by testing \emph{individual} samples with multiple tests, using a covariance scheme\cr
  \code{\link{betaPERT}}: calculate the parameters of a Beta-PERT distribution\cr
  \code{\link{betaExpert}}: calculate the parameters of a Beta distribution based on expert opinion
}

\examples{
## Sandflies in Aurabani, Nepal, 2007
pool_results <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
pool_sizes <- c(2, 1, 6, 10, 1, 7, 1, 4, 1, 3)

## Sensitivity ranges uniformly between 60% and 95%
## Specificity is considered to be 100%

#> BUGS-style:
truePrevPools(x = pool_results, n = pool_sizes,
              SE = ~dunif(0.60, 0.95), SP = 1)
			  
#> list-style:
SE <- list(dist = "uniform", min = 0.60, max = 0.95)
truePrevPools(x = pool_results, n = pool_sizes,
              SE = SE, SP = 1)
}
