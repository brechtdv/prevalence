\name{prev-class}
\Rdversion{1.1}
\docType{class}
\alias{prev-class}

\title{Class \code{"prev"}}
\description{
The \code{"prev"} class represents output from Bayesian true prevalence estimation models.
}

\section{Objects from the Class}{
Objects of class \code{"prev"} are created by \code{\link{truePrev}} and \code{\link{truePrevPools}}.
}

\section{Slots}{
Objects of class \code{"prev"} contain the following four slots:
  \describe{
    \item{\code{par}:}{
	  A list of input parameters
	}
    \item{\code{model}:}{
	  The fitted Bayesian model, in BUGS language (class \code{"prevModel"}) 
	}
    \item{\code{mcmc}:}{
	  A list, with one element per chain, of the simulated true prevalences
	}
    \item{\code{diagnostics}:}{
	  A list with elements for Deviance Information Criterion (\code{$DIC}) and Brooks-Gelman-Rubin diagnostics (\code{$BGR})
	}
  }
}

\author{
\email{Brecht.Devleesschauwer@UGent.be}
}

\seealso{
  \code{\link{truePrev}}, \code{\link{truePrevPools}}\cr
  \code{\link{show-methods}}, \code{\link{print-methods}}, \code{\link{summary-methods}}, \code{\link{plot-methods}}
}

\examples{
## Taenia solium cysticercosis in Nepal
SE <- list(dist = "uniform", min = 0.60, max = 1.00)
SP <- list(dist = "uniform", min = 0.75, max = 1.00)
TP <- truePrev(x = 142, n = 742, SE = SE, SP = SP)

## Summarize estimates per chain
summary(TP)

## Diagnostic plots
par(mfrow = c(2, 2))
plot(TP)

## Use 'str()' to see the structure of object TP
str(TP)

## Every slot can be accessed using the '@' operator
TP@par          # input parameters
TP@model        # fitted model
TP@mcmc         # simulated true prevalences
TP@diagnostics  # DIC and BGR

## TP@mcmc inherits from class 'mcmc.list' in coda package
## List all available methods for this class
methods(class = "mcmc.list")

## See library(help = "coda") for more available functions
## For instance, plot(TP) shows a combination of the following
densplot(TP@mcmc)
traceplot(TP@mcmc)
gelman.plot(TP@mcmc)
autocorr.plot(TP@mcmc)
}

\keyword{classes}
