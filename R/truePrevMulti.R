## ---------------------------------------------------------------------------#
## USER INTERFACE ------------------------------------------------------------#

truePrevMulti <-
function(x, n, prior, method = c("conditional", "covariance"),
         conf.level = 0.95, nchains = 2, burnin = 5000, update = 10000,
         verbose = FALSE) {

  ## check method
  method <- match.arg(method)

  ## check x and n
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  checkInput(x, "x", class = "integer", min = 0)
  checkInput(n, "n", class = "integer", minEq = 0)
  if (sum(x) != n) stop("'x' does not sum to 'n'")
  if ((log(length(x), 2) %% 1 != 0) | length(x) < 4) {
    stop("'x' is not correctly specified; see ?definition_x")
  }

  ## check prior
  if (missing(prior)) stop("'prior' is missing")
  prior <-
    switch(method,
           conditional = checkMultiPrior_conditional(substitute(prior)),
           covariance = checkMultiPrior_covariance(substitute(prior),
                                                   log(length(x), 2)))

  ## check conf.level
  checkInput(conf.level, "conf.level", class = "numeric", range = c(0, 1))

  ## check nchains, burnin & update
  checkInput(nchains, "nchains", class = "integer", min = 2)
  checkInput(burnin, "burnin", class = "integer", min = 1)
  checkInput(update, "update", class = "integer", min = 1)

  ## check options
  checkInput(verbose, "verbose", class = "logical")

  ## get output
  out <-
    switch(method,
           conditional =
             truePrevMultinom_conditional(x, n, prior, conf.level,
                                          nchains, burnin, update, verbose),
           covariance = 
             truePrevMultinom_covariance(x, n, prior, conf.level,
                                         nchains, burnin, update, verbose))

  ## return output
  return(out)
}


## ---------------------------------------------------------------------------#
## CHECK PRIOR FOR CONDITIONAL DEPENDENCE MODEL ------------------------------#

checkMultiPrior_conditional <-
function(prior){
  ## evaluate whether prior is defined correctly
  first_element <- as.character(prior)[1]
  if (!any(c("{", "list") == first_element))
    stop("'prior' is not defined correctly")

  ## if prior is defined as a list
  ## note: list element names currently not taken in account!
  if (first_element == "list"){
    n <- length(prior) - 1
    priors_list <- vector("list", n)
    for (i in seq(n))
      priors_list[[i]] <- checkSeSp(eval(parse(text = prior)[[i + 1]]))
  }

  ## if prior is defined as a function
  if (first_element == "{"){
    n <- length(prior) - 1
    priors_list0 <- vector("list", n)
    for (i in seq(n)) {
      priors_list0[[i]] <- explode(as.character(prior[[i + 1]]), "cond")
    }

    ## get indices from priors_list0
    index <- sapply(priors_list0, function(x) as.numeric(x[[1]]))

    ## check if all indices exist
    if (length(index) != max(index)) {
      stop("The indices of 'theta[.]' are not correctly specified.\n",
           "See ?theta for more info.")
    }
    if (!all(unique(index) == index)) {
      stop("The indices of 'theta[.]' are not correctly specified.\n",
           "See ?theta for more info.")
    }
    if (length(index) == 1 | log(length(index) + 1, 2)%%1 != 0) {
      stop("The number of specified theta values is incorrect.\n",
           "See ?theta for more info.")
    }

    ## re-arrange list elements if needed
    order <- order(index)
    priors_list <- vector("list", n)
    for (i in seq(n))
      priors_list[[i]] <- priors_list0[[order[i]]][[2]]
    
  }

  ## return prior in list format
  return(priors_list)
}


## ---------------------------------------------------------------------------#
## CHECK PRIOR FOR COVARIANCE MODEL ------------------------------------------#

checkMultiPrior_covariance <-
function(prior, h){
  ## evaluate whether prior is defined correctly
  first_element <- as.character(prior)[1]
  if (!any(c("{", "list") == first_element))
    stop("'prior' is not defined correctly")

  ## if prior is defined as a list
  ## note: list element names currently not taken in account!
  if (first_element == "list"){
    n <- length(prior) - 1
    priors_list <- vector("list", n)
    for (i in seq(n))
      priors_list[[i]] <- checkSeSp(eval(parse(text = prior)[[i + 1]]))
  }

  ## if prior is defined as a function
  if (first_element == "{"){
    n <- length(prior) - 1
    priors_list <- vector("list", n)
    for (i in seq(n)) {
      priors_list[[i]] <-
        explode(as.character(prior[[i + 1]]), "covariance")
    }

    ## based on h tests, these priors are expected
    priors <- get_nodes(h)

    ## check if all priors are defined
    priors_nodes <- sapply(priors_list, function(x) x[[1]])
    if (!all(priors_nodes %in% priors)) {
      stop("priors are not correctly specified")
    }

    ## re-arrange list elements if needed
    order <- match(priors_nodes, priors)
    priors_list <- priors_list[order]
  }

  ## return prior in list format
  return(priors_list)
}


## ---------------------------------------------------------------------------#
## MAIN EXPLODE FUNCTION -----------------------------------------------------#

explode <-
function(x, method){
  priors <- vector("list", 2)
  priors[[1]] <-
    switch(method,
           conditional = explode_theta(x[2]),
           covariance = explode_nodes(x[2]))

  explode_operator(x[1])

  type <-
    ifelse(strsplit(priors[[1]], "[", fixed = T)[[1]][1] %in% c("a", "b"),
           "cov", "prob")
  priors[[2]] <- explode_dist(x[3], type)

  return(priors)
}


## ---------------------------------------------------------------------------#
## EXPLODE THETA (CONDITIONAL DEPENDENCE MODEL) ------------------------------#

explode_theta <-
function(x){
  ## find 'theta[]'
  if (length(grep("theta", x, fixed = TRUE)) != 1)
    stop("Priors must be defined as vector 'theta'")
  if (length(grep("[", x, fixed = TRUE)) != 1 |
      length(grep("]", x, fixed = TRUE)) != 1)
    stop("The different values of theta must be defined as 'theta[.]'")

  ## extract '.' in 'theta[.]'
  x <- strsplit(x, "theta[", fixed = TRUE)[[1]][2]
  theta <- strsplit(x, "]", fixed = TRUE)[[1]][1]

  ## theta should be an integer
  if (!is.numeric(theta) && as.numeric(theta) %% 1 != 0)
    stop("'theta[.]' not specified correctly")

  return(theta)
}


## ---------------------------------------------------------------------------#
## EXPLODE NODES (COVARIANCE MODEL) ------------------------------------------#

explode_nodes <-
function (x) {
  if (!any(c(grepl("TP", x, fixed = TRUE),
             grepl("SE", x, fixed = TRUE),
             grepl("SP", x, fixed = TRUE),
             grepl("a", x, fixed = TRUE),
             grepl("b", x, fixed = TRUE)))) {
      stop("Priors must be named 'TP', 'SE', 'SP', 'a' or 'b'")
  }

  if (!(strsplit(x, "[", fixed = T)[[1]][1] %in%
        c("TP", "SE", "SP", "a", "b"))){
      stop("Priors must be named 'TP', 'SE', 'SP', 'a' or 'b'")
  }

  if (x != "TP" &&
      !all(c(grepl("[", x, fixed = TRUE),
             grepl("]", x, fixed = TRUE)))) {
      stop("Priors must be defined as vectors")
  }

  if (x != "TP") {
    rhs <- strsplit(x, "[", fixed = TRUE)[[1]][2]
    i <- strsplit(rhs, "]", fixed = TRUE)[[1]][1]
    if (!grepl("^[[:digit:]]+$", i) || i < 1) {
      stop("Prior '", strsplit(x, "[", fixed = T)[[1]][1],
           "' not correctly indexed")
    }
  }

  return(x)
}


## ---------------------------------------------------------------------------#
## EXPLODE OPERATOR ----------------------------------------------------------#

explode_operator <-
function(operator){
  ## operator should be '~' or '<-' or '='
  if (!any(c("~", "<-", "=") == operator))
    stop("Operator should be either '~', '<-' or '='")
}


## ---------------------------------------------------------------------------#
## EXPLODE DISTRIBUTION ------------------------------------------------------#

explode_dist <-
function(x, type){
  d <- dist2list(x, type)
  return(d)
}


## ---------------------------------------------------------------------------#
## DEFINE NODES OF CONDITIONAL MODEL -----------------------------------------#

get_nodes <-
function(h) {
  nodes <-
    c("TP",
      paste("SE[", seq(h), "]", sep = ""),
      paste("SP[", seq(h), "]", sep = ""),
      paste("a[", seq(sum(choose(h, seq(h, 2)))), "]", sep = ""),
      paste("b[", seq(sum(choose(h, seq(h, 2)))), "]", sep = ""))
  return(nodes)
}


## ---------------------------------------------------------------------------#
## CREATE MODEL FOR METHOD == CONDITIONAL ------------------------------------#

truePrevMultinom_conditional <-
function(x, n, prior, conf.level,
         nchains, burnin, update, verbose) {

  ## create model
  t <- log(length(x), 2)     # number of tests
  ntheta <- 2 ^ (t + 1) - 1  # number of thetas

  model <- character()

  ## write model initiation
  model[1] <- "model {"
  model[2] <- paste("x[1:", 2 ^ t,
                    "] ~ dmulti(AP[1:", 2 ^ t, "], n)",
                    sep = "")

  ## write AP[] definitions in terms of theta[]
  s <- multiModel_select(t)  # define theta construct for SE/SP
  p <- multiModel_probs(s)   # define AP[.] in terms of theta[.]
  model <- c(model, "", p, "")

  ## write theta[] prior
  for (i in seq(ntheta))
    model <- c(model,
      writeSeSp(paste("theta[", i, "]", sep = ""), prior[[i]]))

  ## write bayesP definition
  bayesP <-
    c(paste("x2[1:", (2^t), "] ~ dmulti(AP[1:", (2^t), "], n)", sep=""),
      paste("for (i in 1:", (2^t), ") {", sep=""),
      "d1[i] <- x[i] * log(max(x[i],1) / (AP[i]*n))",
      "d2[i] <- x2[i] * log(max(x2[i],1) / (AP[i]*n))",
      "}",
      "G0 <- 2 * sum(d1[])",
      "Gt <- 2 * sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  model <- c(model, "", bayesP)

  ## write SE[]/SP[] definition
  model <- c(model, "", multiModel_SeSp(t))

  ## close model
  model <- c(model, "}")

  ## define model class
  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, n = n)

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  nodes <- paste(c("SE", "SP"), rep(seq(t), each = 2), sep = "")
  nodes <- c("TP", nodes, "bayesP")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = nodes, verbose = verbose)

  ## define mcmc samples
  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")
  names <- colnames(mcmc.list[[1]])
  mcmc.list_list <- list()
  order <- c(length(names) - 1, c(t(cbind(1:t, 1:t+t))), length(names))
  for (i in seq_along(names))
    mcmc.list_list[[i]] <- mcmc.list[, order[i]]
  names(mcmc.list_list) <- names[order]

  ## define diagnostics
  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list_list$TP, autoburnin = FALSE)$psrf)
  bayesP <- mean(unlist(mcmc.list_list$bayesP))

  ## define output
  out <-
    new("prev",
        par = list(x = x, n = n, prior = prior,
                   conf.level = conf.level, nchains = nchains,
                   burnin = burnin, update = update, inits = inits),
        model = model,
        mcmc = mcmc.list_list,
        diagnostics = list(DIC = DIC,
                           BGR = data.frame(mean = BGR[1],
                                            upperCL = BGR[2]),
                           bayesP = bayesP))

  ## return output
  return(out)
}


## ---------------------------------------------------------------------------#
## CREATE MODEL FOR METHOD == COVARIANCES ------------------------------------#

truePrevMultinom_covariance <-
function(x, n, prior, conf.level,
         nchains, burnin, update, verbose) {

  ## number of tests
  h <- log(length(x), 2)

  ## number of priors
  n_priors <- 1 + (2 * h) + (2 * sum(choose(h, seq(h, 2))))

  ## define model vector
  model <- character()

  ## write model initiation
  model[1] <- "model {"
  model[2] <- paste("x[1:", 2 ^ h,
                    "] ~ dmulti(AP[1:", 2 ^ h, "], n)",
                    sep = "")

  ## write prob_se[], prob_sp[]
  s <- multiModel_select(h)[[1]]
  s <- s[rev(seq(nrow(s))), ]

  prob_se <- paste("prob_se[", seq(nrow(s)), "] <-", sep = "")

  for (i in seq(nrow(s))) {
    ## first element
    prob_se[i] <-
      paste(prob_se[i],
            paste(ifelse(s[i, ] == 1,
                         paste("SE[", seq(ncol(s)), "]", sep = ""),
                         paste("(1 - SE[", seq(ncol(s)), "])", sep = "")),
                  collapse = " * "))

    ## define index for 'a'
    a <- c(1, 0)
    
    ## h - 2 elements
    if (h > 2) {
      for (k in seq((h - 2), 1)) {
        se <-
          apply(t(apply(combn(h, k), 1, rev)),
                2,
                function(x) {
                  paste(ifelse(s[i, x] == 1,
                               paste("SE[", x, "]", sep = ""),
                               paste("(1 - SE[", x, "])", sep = "")),
                        collapse = " * ")
                })

        sign <-
          apply(t(apply(combn(h, k), 1, rev)),
                2,
                function(x) prod(2 * s[i, -x] - 1))  # convert (1,0) to (1,-1)

        a[2] <- a[2] + choose(h, k)
        prob_se[i] <-
          paste(prob_se[i],
                paste(
                  paste(ifelse(sign == 1,
                               " + ", " - "),
                        "a[", seq(a[1], a[2]), "] * ", se,
                        sep = ""),
                  collapse = ""),
                 sep = "")
        a[1] <- a[2] + 1
      }
    }

    ## final element
    prob_se[i] <-
      paste(prob_se[i],
            paste(ifelse(prod(2 * s[i, ] - 1) == 1,
                  "+ ","- "),
                  "a[", a[1], "]",
                  sep = ""))
  }

  prob_sp <- gsub("SE", "SP", prob_se)
  prob_sp <- gsub("prob_se", "prob_sp", prob_sp)
  prob_sp <- gsub("a", "b", prob_sp)
  for (i in seq_along(prob_sp)) {
    prob_sp[i] <-
      gsub(paste("prob_sp[", i ,"]", sep = ""),
           paste("prob_sp[", 1 + length(prob_sp) - i ,"]", sep = ""),
           fixed = TRUE,
           prob_sp[i])
  }

  model <-
    c(model,
      "",
      prob_se,
      "",
      prob_sp, "")

  ## write definition of AP and constraints
  model <-
    c(model,
      paste("for (i in 1:", 2 ^ h, ") {", sep = ""),
      "AP[i] <- TP * prob_se[i] + (1 - TP) * prob_sp[i]",
      "",
      write_constraint("AP", 1, 1),
      write_constraint("AP", 2, 2),
      write_constraint("prob_se", 1, 3),
      write_constraint("prob_se", 2, 4),
      write_constraint("prob_sp", 1, 5),
      write_constraint("prob_sp", 2, 6),
      "}",
      "")

  ## write prior
  priors <- get_nodes(h)

  for (i in seq(n_priors)) {
    model <-
      c(model,
        writeSeSp(prior[[i]][[1]], prior[[i]][[2]]))
  }

  ## write Bayes-P definition
  model <- c(model, "", write_bayesP(h))

  ## close model
  model <- c(model, "}")

  ## define model class
  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, x2 = x, n = n,
               O1 = rep(1, 2 ^ h), O2 = rep(0, 2 ^ h),
               O3 = rep(1, 2 ^ h), O4 = rep(0, 2 ^ h),
               O5 = rep(1, 2 ^ h), O6 = rep(0, 2 ^ h))

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  nodes <- c("TP", "SE", "SP", "a", "b", "bayesP")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = nodes, verbose = verbose)

  ## define mcmc samples
  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")
  names <- colnames(mcmc.list[[1]])
  if (h == 2) {
    names[which(names == "a")] <- "a[1]"
    names[which(names == "b")] <- "b[1]"
  }
  mcmc.list_list <- list()
  order <- match(c(priors, "bayesP"), names)
  for (i in seq_along(names))
    mcmc.list_list[[i]] <- mcmc.list[, order[i]]
  names(mcmc.list_list) <- names[order]

  ## define diagnostics
  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list_list$TP, autoburnin = FALSE)$psrf)
  bayesP <- mean(unlist(mcmc.list_list$bayesP))

  ## define output
  out <-
    new("prev",
        par = list(x = x, n = n, prior = prior,
                   conf.level = conf.level, nchains = nchains,
                   burnin = burnin, update = update, inits = inits),
        model = model,
        mcmc = mcmc.list_list,
        diagnostics = list(DIC = DIC,
                           BGR = data.frame(mean = BGR[1],
                                            upperCL = BGR[2]),
                           bayesP = bayesP))

  ## return output
  return(out)
}


## ---------------------------------------------------------------------------#
## WRITE DEFINITION OF BAYES-P -----------------------------------------------#

write_bayesP <-
function(h) {
  bayesP <-
    c(paste("x2[1:", (2^h), "] ~ dmulti(AP[1:", (2^h), "], n)", sep = ""),
      paste("for (i in 1:", (2^h), ") {", sep = ""),
      "d1[i] <- pow(x[i] - AP[i] * n, 2) / (AP[i] * n)",
      "d2[i] <- pow(x2[i] - AP[i] * n, 2) / (AP[i] * n)",
      "}",
      "G0 <- sum(d1[])",
      "Gt <- sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  return(bayesP)
}

write_bayesP0 <-
function(h) {
  bayesP <-
    c(paste("x2[1:", (2^h), "] ~ dmulti(AP[1:", (2^h), "], n)", sep = ""),
      paste("for (i in 1:", (2^h), ") {", sep = ""),
      "d1[i] <- x[i] * log(max(x[i],1) / (AP[i]*n))",
      "d2[i] <- x2[i] * log(max(x2[i],1) / (AP[i]*n))",
      "}",
      "G0 <- 2 * sum(d1[])",
      "Gt <- 2 * sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  return(bayesP)
}


## ---------------------------------------------------------------------------#
## WRITE CONSTRAINT ON PROB_SE, PROB_SP --------------------------------------#

write_constraint <-
function(node, constraint, i) {
  add <- ifelse (constraint == 2, " - 1", "")
  constr <-
    c(paste("constraint", i, "[i] <- step(", node, "[i]", add, ")", sep = ""),
      paste("O", i, "[i] ~ dbern(constraint", i, "[i])", sep = ""))
  return(constr)
}