truePrevBinom <-
function(x, n, Se, Sp, prior,
         nchains, burnin, update, verbose){

  ## create model
  model <- character()
  model[1] <- "model {"
  model[2] <- "x ~ dbin(AP, n)"
  model[3] <- "AP <- SE * TP + (1 - SP) * (1 - TP)"

  model <- c(model, writeSeSp("SE", Se))
  model <- c(model, writeSeSp("SP", Sp))

  model <-
    c(model,
      paste0("TP ~ dbeta(", prior[1], ", ", prior[2], ")"))

  model <- c(model, "}")

  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, n = n)

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = c("TP", "SE", "SP"), verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  ## deviance information criterion
  DIC <- JAGSout$dic

  ## brooks-gelman-rubin diagnostic
  is_stoch <- apply(mcmc.list[[1]], 2, var) > 0
  mcmc.list_stoch <- mcmc.list[, is_stoch]
  BGR <- c(gelman.diag(mcmc.list_stoch, autoburnin = FALSE)$psrf)

  ## create new 'prev' object
  out <-
    new("prev",
        par = list(x = x, n = n, SE = Se, SP = Sp, prior = prior,
                   nchains = nchains, burnin = burnin, update = update,
                   inits = inits),
        model = model,
        mcmc = mcmc.list,
        diagnostics = list(DIC = DIC,
                           BGR = data.frame(mean = BGR[1], upperCL = BGR[2])))

  return(out)
}