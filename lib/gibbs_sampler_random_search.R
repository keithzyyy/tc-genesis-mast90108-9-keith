set.seed(123)

### Libraries ----
if (!require(tidyr)) install.packages("tidyr")
if (!require(coda)) install.packages("coda")
if (!require(zoo)) install.packages("zoo")
if (!require(cowplot)) install.packages("cowplot")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(coda)    # effectiveSize
  library(zoo)     # rollmean for smoothing
  library(cowplot) # for grid plotting
})


#' computes AIC of the logistic regression model comprised
#' of predictors indexed by indicator variable v
#'
#' @param data = the data frame comprised of cbind(X,y) where
#' X(nxp) is the data matrix and y(nx1) is the response. 
#' 
#' @param v = the current state in the markov chain.
#' 
#' @param method = either use glm() or manually compute the AIC.
#' 
#' @returns AIC of the model indexed by v 
compute.AIC = function(data,v,debug=FALSE){
  
  # there may be some components of v that are 0,
  # meaning we should exclude them from the model.
  # note data is (n x (p)) due to response being in 
  # the last column
  X.sub = data[, which(v == 1), drop=FALSE]
  y.sub = data[, ncol(data)]
  
  #print(v)
  
  data.model = data.frame(X.sub, y = y.sub)
  
  # Check for degenerate case
  # if (ncol(X.sub) == 0 || length(unique(y.sub)) == 1) return(Inf)
  
  # intercept is NOT included in the variable selection -- 
  # so is always included in the model.
  #fit = tryCatch({
  #  mod = glm(y ~ ., data = data.model, family = binomial)
  #  #DIAGNOSTICS
  #  if (debug){
  #    #cat("Dimension of data model: ", dim(data.model),"\n")
  #    #cat("Number of parameters: ", length(mod$coefficients), "\n")
  #    #print(summary(mod))
  #  }
  #}, error = function(e) return(NULL))
  
  #if (is.null(fit) || is.null(fit$coefficients) || any(is.na(fit$coefficients))) {
  #  print("WARNING: Possible Singularities, AIC returns Inf")
  #  return(Inf)
  #}
  
  # fit a glm model 
  fit = glm(y ~ ., data = data.model, family = binomial)
  
  # Compute probabilities
  X.sub.intercept = cbind(rep(1,nrow(data.model)), X.sub)
  eta = as.vector(as.matrix(X.sub.intercept) %*% fit$coefficients)
  probs = 1 / (1 + exp(-eta))
  
  # Clip probs to avoid log(0)
  eps = 1e-10
  probs = pmin(pmax(probs, eps), 1 - eps)
  
  # Compute NLL
  nll = -sum(y.sub * log(probs) + (1 - y.sub) * log(1 - probs))
  
  # AIC = 2 * NLL + 2 * d
  d = sum(v == 1)
  AIC.mod = nll + d
  
  # cat("AIC for model with covariate idxs ", which(v==1), " is ", AIC.mod, "\n")
  return(AIC.mod)
  
}



#' computes BIC of the logistic regression model comprised
#' of predictors indexed by indicator variable v
#'
#' @param data = the data frame comprised of cbind(X,y) where
#' X(nxp) is the data matrix and y(nx1) is the response
#' 
#' @param v = the current state in the markov chain.
#' 
#' @returns BIC of the model indexed by v
compute.BIC = function(data,v){
  
  AIC.mod = compute.AIC(data,v)
  #from the HK paper,
  #AIC = NLL + p, and BIC = NLL + 0.5*p*log(N)
  N = nrow(data)
  p = sum(v==1)
  
  return(AIC.mod - p + 0.5*p*log(N))
}

#' Simulates an indicator RV from the conditional distribution of
#' the indicator variables.
#'
#'
#' @param k = the column index. Must be from 1 to p
#' @param data = a data frame comprised of a X(nxp) design matrix
#' and a y(nx1) response vector. Should be in the format: cbind(X,y)
#' @param v = the current state in the markov chain.
#' @param penalty = penalty type of interest to compute the model criterion
#' 
#' @returns simulated value from the conditional distribution (binary)
#simulate.conditional.vk = function(k, data, v,
#                                   penalty,
#                                   debug=TRUE,
#                                   lambda=1){

# Want: P(v) / [P(v0) + P(v1)]
# v0 := v with v_k = 0
# v1 := v with v_k = 1
#  v0 <- v; v0[k] <- 0
#  v1 <- v; v1[k] <- 1


# obtain the success probability P(v_k=1|v-k), given by:
# p.succ = prob.joint / (prob.joint.0 + prob.joint.1)
#       = exp(-criterion(v1)) / [ exp(-criterion(v0)) + exp(-criterion(v1)) ] 

#  if (penalty == 'AIC'){
#    log.p0 = -compute.AIC(data, v0) * lambda
#    log.p1 = -compute.AIC(data, v1) * lambda
#  }

#  if (penalty == 'BIC'){
#    log.p0 = -compute.BIC.auto(data, v0) * lambda
#    log.p1 = -compute.BIC.auto(data, v1) * lambda
#  }


# Use log-sum-exp trick to avoid underflow
#  max.log.p = max(log.p0, log.p1)
#  denom = log(exp(log.p0 - max.log.p) + exp(log.p1 - max.log.p)) + max.log.p
#  log.p.succ = log.p1 - denom
#  p.succ = exp(log.p.succ)
#  v.sim = rbinom(1,1,p.succ)

#  if (debug){

#cat("\nv0=",v0,"\n")

#    cat("Updating coordinate ", k, "\n")
#    cat("p.succ=",p.succ,"\n")
#    cat("simulated value=",v.sim,"\n")

#cat("\nlog.p.succ=",log.p.succ,"\n")

#  }

#  if (p.succ == 0 || is.na(p.succ)){
#    return(0)
#  }

# finally simulate from associated bernoulli dist
#  return(v.sim)

#}

simulate.conditional.vk <- function(k, data, v,
                                    penalty,
                                    lambda = 1,
                                    debug  = FALSE) {
  v0 <- v; v0[k] <- 0
  v1 <- v; v1[k] <- 1
  
  B0 <- if (penalty == "AIC") compute.AIC(data, v0) else compute.BIC.auto(data, v0)
  B1 <- if (penalty == "AIC") compute.AIC(data, v1) else compute.BIC.auto(data, v1)
  
  # p(v_k=1 | v_-k) = 1 / (1 + exp(lambda*(B1 - B0)))
  p1 <- plogis(-lambda * (B1 - B0))   # numerically stable
  
  # edge case: if incl prob is extreme return 1/0
  if (is.na(p1) || p1 <= 0){
    vk = 0
  } else if (p1 >= 1){
    vk = 1
  } else {
    # otherwise sample from cond probabilility dist
    vk = rbinom(1L, 1L, p1)
  }
  
  if (debug) {
    cat(sprintf("k=%d  B0=%.3f  B1=%.3f  p1=%.6f, v_k=%d\n", k, B0, B1, p1, vk))
  }
  
  vk
}




compute.BIC.auto <- function(data, v) {
  y <- data[, ncol(data)]
  if (!any(v == 1)) {
    # Intercept-only logistic model
    fit <- glm(y ~ 1, data = data.frame(y = y), family = binomial())
  } else {
    X.sub <- data[, which(v == 1), drop = FALSE]
    fit   <- glm(y ~ ., data = data.frame(X.sub, y = y), family = binomial())
  }
  BIC(fit)
}




MH.random.search <- function(v0, data,
                             burnin = 200,
                             iter   = 1200,
                             debug  = TRUE,
                             lambda = 1) {
  stopifnot(iter > 1, burnin >= 0, iter > burnin)
  start.time <- Sys.time()
  
  if (burnin == 0) cat("burnin is 0 → continuing an existing chain.\n")
  
  p        <- ncol(data) - 1
  samples  <- matrix(0, iter, p)
  samples[1, ] <- v0
  acc      <- 0L
  
  # init BIC at v0 (cached)
  bic_curr <- compute.BIC.auto(data, v0)
  
  # (optional) diagnostics
  bic_trace <- numeric(iter); bic_trace[1] <- bic_curr
  best_bic  <- bic_curr; best_v <- v0
  
  if (!requireNamespace("progress", quietly = TRUE)) install.packages("progress")
  pb <- progress::progress_bar$new(format = " Iter [:bar] :percent eta: :eta",
                                   total = iter, clear = FALSE, width = 60)
  
  for (j in 2:iter) {
    v.prev <- samples[j - 1, ]
    v      <- v.prev
    
    # proposal = flip one bit
    q <- sample.int(p, 1L)       
    v[q] <- 1L - v.prev[q]
    
    bic_prop <- compute.BIC.auto(data, v)
    dBIC     <- bic_prev_minus_prop <- bic_curr - bic_prop  # = -(BIC_prop - BIC_curr)
    
    accepted <- (log(runif(1)) <= lambda * dBIC)
    if (accepted) {
      samples[j, ] <- v
      bic_curr     <- bic_prop
      acc          <- acc + 1L
    } else {
      samples[j, ] <- v.prev
    }
    
    # track best & trace
    bic_trace[j] <- bic_curr
    if (bic_curr < best_bic) { best_bic <- bic_curr; best_v <- samples[j, ] }
    
    
    
    if (debug) {
      rel_freq <- mean(samples[1:j, q])
      status   <- if (accepted) "ACCEPTED" else "REJECTED"
      cat(sprintf("Iter %d | Var=%s | %s | rel.freq=%.2f | #active=%d | acc.rate=%.3f | bic value=%.3f\n",
                  j,
                  colnames(data)[q],
                  status,
                  rel_freq,
                  sum(samples[j, ]),
                  acc / j,
                  bic_curr)
      )
    }
    
    pb$tick()
  }
  
  colnames(samples) <- colnames(data)[1:p]
  
  samples.post <- if (burnin > 0) samples[-seq_len(burnin), , drop = FALSE] else samples
  
  end.time <- Sys.time()
  cat("\nTime taken (min):", as.numeric(end.time - start.time, units = "mins"), "\n")
  
  list(
    samples   = samples.post,
    acc.rate  = acc / (iter - 1),
    bic_trace = bic_trace,
    best = list(v = best_v, BIC = best_bic)
  )
}




Gibbs.random.search.lambda <- function(v0, data,
                                       burnin  = 200,
                                       iter    = 800,
                                       penalty = "BIC",
                                       debug   = FALSE,
                                       lambda  = 1) {
  stopifnot(iter > 1, burnin >= 0, iter > burnin)
  
  start.time <- Sys.time()
  
  p <- ncol(data) - 1L
  xnames <- colnames(data)[1:p]
  
  # storage
  samples   <- matrix(0L, nrow = iter, ncol = p)
  samples[1, ] <- v0
  bic_trace <- numeric(iter)
  
  # initial BIC
  bic_curr <- if (penalty == "AIC") compute.AIC(data, v0) else compute.BIC.auto(data, v0)
  bic_trace[1] <- bic_curr
  best_bic <- bic_curr
  best_v   <- v0
  best_it  <- 1L
  
  # optional progress
  if (debug) {
    if (!requireNamespace("progress", quietly = TRUE)) install.packages("progress")
    pb <- progress::progress_bar$new(
      format = " Iter [:bar] :percent eta: :eta",
      total = iter, clear = FALSE, width = 60
    )
  }
  
  for (j in 2:iter) {
    cat(sprintf("\n [%s]  Begin iteration %d out of %d \n", Sys.time(), j, iter))
    v <- samples[j - 1L, ]
    
    # one full Gibbs sweep over coordinates
    for (k in 1:p) {
      v[k] <- simulate.conditional.vk(k, data, v,
                                      penalty = penalty,
                                      lambda  = lambda,
                                      debug   = debug)
    }
    
    samples[j, ] <- v
    
    if (debug){
      cat("Number of variables in v: ", sum(v), "\n")
    }
    
    bic_curr <- if (penalty == "AIC") compute.AIC(data, v) else compute.BIC.auto(data, v)
    bic_trace[j] <- bic_curr
    
    if (bic_curr < best_bic) {
      best_bic <- bic_curr
      best_v   <- v
      best_it  <- j
    }
    
    if (debug) pb$tick()
  }
  
  colnames(samples) <- xnames
  
  # burn-in handling
  post_samples <- if (burnin > 0) samples[-seq_len(burnin), , drop = FALSE] else samples
  post_trace   <- if (burnin > 0) bic_trace[-seq_len(burnin)] else bic_trace
  
  # helpers for best model
  best_key  <- paste0(best_v, collapse = "")
  best_vars <- xnames[as.logical(best_v)]
  
  cat("\nTime taken (min):", round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2), "\n")
  
  list(
    samples   = post_samples,
    bic_trace = post_trace,
    best = list(
      v        = best_v,
      BIC      = best_bic,
      iter     = best_it,
      key      = best_key,
      variables = best_vars
    )
  )
}





#' Performs the vanilla MCMC random search procedure to find the model
#' that minimizes some model selection criterion.
#'
#'
#' @param data = a matrix/dataframe which is a concatenation of 2 things:
#' - X = the n x p design matrix
#' - y = the response vector (a nx1 matrix)
#' @param v0 = initial values for p-dimensional indicator vectors
#'
#' @returns a matrix of size (iter-burnin) x (p), where p is the 
#' dimension of beta 

Gibbs.random.search = function(v0, data,
                               burnin=200,
                               iter=800,
                               penalty='AIC',
                               debug=TRUE){
  
  start.time.function = Sys.time()
  
  if (burnin == 0){
    cat("Burnin is set to 0. CONTINUING AN EXISTING CHAIN.\n")
  }
  
  # 0. defining essential quantities
  p = ncol(data) - 1
  samples = matrix(0, iter, p)
  samples[1,] = v0
  
  # progress bar
  if (!require(progress)) install.packages("progress")
  library(progress)
  pb <- progress_bar$new(
    format = "  Iteration [:bar] :percent eta: :eta",
    total = iter, clear = FALSE, width = 60
  )
  
  
  # 1. begin the gibbs sampler iteration!
  for (j in 2:iter){ # j = MCMC iteration idx
    
    # 1.0 logging
    #if ((j-2) %% 50 == 0){
    #  cat(sprintf("\n [%s]  Begin iteration %d \n", Sys.time(), j))
    #}
    cat(sprintf("\n [%s]  Begin iteration %d out of %d \n", Sys.time(), j, iter))
    
    # 1.1 obtain the most recent state of the markov chain
    # v = (v^(j-1)_1, v^(j-1)_2, ... , v^(j-1)_p)
    v.prev = samples[j-1,]
    
    # initialize the next state, which will be updated component wise
    v = v.prev
    
    if (debug){
      cat("Number of variables in v.prev: ", sum(v.prev), "\n")
    }
    
    # 1.2 simulate next state via univariate conditional distributions
    for (k in 1:p){ # k = covariate idx
      # simulate v^(j)_k from P(v_k | all v^(j-1)_{-k})
      v.k =  simulate.conditional.vk(k, data, v, penalty)
      
      # update k-th component of v from 
      # v^(j-1)_{k} to v^(j)_{k}
      v[k] = v.k
    }
    
    # 1.3 we can store the next state. 
    samples[j,] = v
    
    # update progress bar.
    pb$tick()
    
  }
  
  
  covariate.names = colnames(data)[1:p]  # assuming last column is response
  colnames(samples) = covariate.names
  
  end.time.function = Sys.time()
  cat("\n")
  time.diff = end.time.function - start.time.function
  cat("\n Overall time taken (minutes): ",
      as.numeric(time.diff, units = "mins"),
      "\n")
  
  # 2. remove first few samples as burn in, and return the result
  if (burnin > 0){
    return(samples[-c(1:burnin), ])
  }
  
  return(samples)
  
}




Gibbs.random.search.RCP = function(v0, data,
                                   burnin=200,
                                   iter=800,
                                   penalty='AIC'){
  
  start.time.function = Sys.time()
  
  # 0. defining essential quantities
  p = ncol(data) - 1
  samples = matrix(0, iter, p)
  samples[1,] = v0
  
  # progress bar
  if (!require(progress)) install.packages("progress")
  library(progress)
  pb <- progress_bar$new(
    format = "  Iteration [:bar] :percent eta: :eta",
    total = iter, clear = FALSE, width = 60
  )
  
  
  # 1. begin the gibbs sampler iteration!
  for (j in 2:iter){ # j = MCMC iteration idx
    
    # 1.0 logging
    if (j %% 50 == 0){
      cat(sprintf("\n [%s]  Begin iteration %d \n", Sys.time(), j))
    }
    
    # 1.1 obtain the most recent state of the markov chain
    # v = (v^(j-1)_1, v^(j-1)_2, ... , v^(j-1)_p)
    v.prev = samples[j-1,]
    
    # initialize the next state, which will be updated component wise
    v = v.prev
    print(v)
    
    # 1.2 simulate next state via univariate conditional distributions
    for (k in 1:p){ # k = covariate idx
      # simulate v^(j)_k from P(v_k | all v^(j-1)_{-k})
      v.k =  simulate.conditional.vk(k, data, v, penalty)
      
      # update k-th component of v from 
      # v^(j-1)_{k} to v^(j)_{k}
      v[k] = v.k
    }
    
    # 1.3 we can store the next state. 
    samples[j,] = v
    
    # update progress bar.
    pb$tick()
    
  }
  
  
  covariate.names = colnames(data)[1:p]  # assuming last column is response
  colnames(samples) = covariate.names
  
  end.time.function = Sys.time()
  cat("\n")
  time.diff = end.time.function - start.time.function
  cat("\n Overall time taken (minutes): ",
      as.numeric(time.diff, units = "mins"),
      "\n")
  
  # 2. remove first few samples as burn in, and return the result
  return(samples[-c(1:burnin), ])
}



#' Performs the MCMC random search procedure to find the model
#' that minimizes some model selection criterion on tropical cyclone
#' dataset, with `basin` variables controlled (not partake in variable
#' selection) by constraining its indicator to be 1
#'
#'
#' @param data = a matrix which is a concatenation of 2 things:
#' - X = the n x p design matrix
#' - y = the response vector (a nx1 matrix)
#' @param v0 = initial values for p-dimensional indicator vectors, must be
#' on the SAME order as X columns in 'data'
#'
#' @returns a matrix of size (iter-burnin) x (p), where p is the 
#' dimension of beta

Gibbs.random.search.basin = function(v0, data,
                                     burnin=200,
                                     iter=800,
                                     penalty='AIC'){
  
  start.time.function = Sys.time()
  
  
  # IMPORTANT: constrain basin indicator variables
  # to be 1.
  basin.idx <- which(grepl("basin", colnames(data)))
  if (!all(v0[basin.idx] == 1)) {
    warning("Basin indicators in v0 were not all set to 1. Forcing them to 1.")
    v0[basin.idx] <- 1
  }
  
  
  # 0. defining essential quantities
  p = ncol(data) - 1
  samples = matrix(0, iter, p)
  samples[1,] = v0
  
  
  # 1. begin the gibbs sampler iteration!
  for (j in 2:iter){ # j = MCMC iteration idx
    
    # 0. logging
    if (j %% 50 == 0){
      cat("Begin iteration ", j, "\n")
    }
    
    # 1.1 obtain the most recent state of the markov chain
    # v = (v^(j-1)_1, v^(j-1)_2, ... , v^(j-1)_p)
    v.prev = samples[j-1,]
    
    # initialize the next state, which will be updated component wise
    v = v.prev
    
    # DIAGNOSTICS
    #ITER_CUTOFF=50
    #if (j <= ITER_CUTOFF){
    #  print(v)
    #}
    
    # 1.2 simulate next state via univariate conditional distributions
    for (k in setdiff(1:p, basin.idx)){ # we only update non-basin variables
      # simulate v^(j)_k from P(v_k | all v^(j-1)_{-k})
      v.k =  simulate.conditional.vk(k, data, v, penalty)
      
      # update k-th component of v from 
      # v^(j-1)_{k} to v^(j)_{k}
      v[k] = v.k
    }
    
    # 1.3 we can store the next state. 
    samples[j,] = v
    
    
    # DIAGNOSTICS
    #if (j <= ITER_CUTOFF){
    #  return(NAN)
    #}
    
    
  }
  
  covariate.names = colnames(data)[1:p]  # assuming last column is response
  colnames(samples) = covariate.names
  
  end.time.function = Sys.time()
  cat("\n")
  time.diff = end.time.function - start.time.function
  cat("Overall time taken (minutes): ", as.numeric(time.diff, units = "mins"))
  
  
  # 2. remove first few samples as burn in, and return the result
  return(samples[-c(1:burnin), ])
  
  
}



#' Given output from Gibbs.random.search.basin(), obtain relative frequencies
#' of each indicator variable appearing in the MCMC sample
#' 
#'
#'
#' @param samples = a matrix of indicator samples stored row-by-row, with
#' a 'basin.idx' attribute. MUST be the return value of Gibbs.random.search.basin()
#' 
#' @param CUTOFF = threhold for which to pick variables from
#'
#' @returns a list of 2 things: relative frequencies of attributes in 
rel.freq.inference = function(samples, CUTOFF=0.5){
  
  # Get all variable names
  var.names <- colnames(samples)
  
  # Identify basin-related variables
  is.basin <- grepl("basin", var.names)
  
  # compute relative frequencies of each indicator variable in the model
  rel.freqs <- colMeans(samples)
  
  # exclude basin variables from analysis
  rel.freqs.non.basin <- rel.freqs[!is.basin]
  
  # select variables above threshold
  selected.vars <- names(rel.freqs.non.basin)[rel.freqs.non.basin >= CUTOFF]
  
  # return relative freqs of indicator variables,
  # as well as the selected variables.
  return(list(
    rel.freqs.non.basin = rel.freqs.non.basin,
    selected.variables = selected.vars
  ))
  
}


# should accept return value from Gibbs.random.search()
find_MAP = function(indicator.samples){
  
  # Convert each row into a string to count identical rows
  # e.g. 1,0,1,0 to 1010
  model.strings <- apply(indicator.samples, 1, paste, collapse = "")
  model.table <- table(model.strings)
  
  # Get the most frequent model (MAP estimate)
  MAP.model.string <- names(which.max(model.table))
  
  # Convert back to binary vector
  MAP.model <- as.numeric(strsplit(MAP.model.string, "")[[1]])
  
  nunique_models = length(model.table)
  
  funique_models = length(model.table) / length(indicator.samples)
  
  cat("\n Number of unique models visited: ",
      nunique_models, "\n")
  cat("\n Proportion of unique models visited out of all draws: ",
      round(100*funique_models, 2), "%\n")
  
  return(MAP.model)
  
}


# Helper: stringify a 0/1 vector to use as a model key
.model_key <- function(v) paste(which(v == 1), collapse = ",")

#' Find unique BIC minimizers (ties kept)
#' @param data  data.frame or matrix cbind(X, y)  (y is last column)
#' @param indicator.samples  K x p binary matrix of inclusion vectors
#' @param top_n  number of *unique BIC levels* to return (ties kept)
#' @param tol  numeric tolerance for declaring BIC ties
#' @return list(models_by_bic, summary, details)
find_BIC_minimizer <- function(data, indicator.samples,
                               top_n = 5, tol = 1e-8,
                               auto_compute_BIC=TRUE) {
  start_time <- Sys.time()
  p <- ncol(indicator.samples)
  xnames <- colnames(data)[seq_len(ncol(data) - 1L)]
  if (is.null(xnames)) xnames <- paste0("X", seq_len(p))
  
  # Compute BIC for every sampled model + build keys
  if (auto_compute_BIC){
    BIC.values <- apply(indicator.samples, 1, function(v) compute.BIC.auto(data, v))
  } else {
    BIC.values <- apply(indicator.samples, 1, function(v) compute.BIC(data, v))
  }
  
  
  
  keys <- apply(indicator.samples, 1, .model_key)
  
  # Collapse duplicate models (same key) — keep the first BIC (it’s deterministic)
  df <- data.frame(
    key   = keys,
    BIC   = BIC.values,
    rowid = seq_along(BIC.values),
    stringsAsFactors = FALSE
  )
  df <- df[!duplicated(df$key), , drop = FALSE]
  
  # Reconstruct variable names per unique model
  uniq_V <- indicator.samples[df$rowid, , drop = FALSE]
  vars_list <- apply(uniq_V, 1, function(v) xnames[as.logical(v)])
  
  # Sort by BIC ascending
  ord <- order(df$BIC)
  df  <- df[ord, , drop = FALSE]
  vars_list <- vars_list[ord]
  
  # Group by BIC within tolerance
  grp <- integer(nrow(df))
  g <- 1L
  grp[1] <- g
  for (i in 2:nrow(df)) {
    if (abs(df$BIC[i] - df$BIC[i - 1]) <= tol) {
      grp[i] <- g
    } else {
      g <- g + 1L
      grp[i] <- g
    }
  }
  df$group <- grp
  
  # Keep only the first top_n unique BIC groups
  keep_groups <- unique(df$group)[seq_len(min(top_n, length(unique(df$group))))]
  keep_idx <- df$group %in% keep_groups
  df_keep <- df[keep_idx, , drop = FALSE]
  vars_keep <- vars_list[keep_idx]
  
  # Build mapping: BIC group -> all tied models (as variable-name vectors)
  models_by_bic <- lapply(split(seq_len(nrow(df_keep)), df_keep$group), function(idx) {
    list(
      BIC   = mean(df_keep$BIC[idx]),  # same within tol
      models = lapply(idx, function(i) vars_keep[[i]])
    )
  })
  
  # Pretty summary table
  summary_tbl <- tibble::tibble(
    Rank      = match(df_keep$group, keep_groups),
    BIC       = df_keep$BIC,
    Variables = vapply(vars_keep, function(v) paste(v, collapse = ", "), character(1)),
    ModelsInThisBIC = ave(df_keep$group, df_keep$group, FUN = length)
  )
  summary_tbl <- summary_tbl[order(summary_tbl$Rank, summary_tbl$BIC), ]
  
  end_time <- Sys.time()
  cat("\nUnique models:", nrow(df), 
      "| Unique BIC levels:", length(unique(df$group)),
      "| Returned levels:", length(keep_groups),
      "| Time (mins):", as.numeric(difftime(end_time, start_time, units = "mins")), "\n")
  
  list(
    models_by_bic = models_by_bic,   # list per BIC level, each with all tied models
    summary = summary_tbl,           # tidy rows; repeated rows for ties within a level
    details = list(all_unique = df, xnames = xnames, tol = tol)
  )
}



#’ Plot running inclusion‐probability (“p‐chart”) for selected predictors
#’
#’ Given a matrix of 0/1 inclusion indicators from an MCMC random‐search chain,
#’ this function plots the cumulative inclusion probability over iterations
#’ for a specified subset of variables (by default, X1–X4).
#’
#’ @param indicator.samples Numeric matrix of size (iterations × p), where each
#’   entry is 0/1 indicating exclusion/inclusion of that predictor in the chain.
#’   Column names should correspond to variable names (e.g. “X1”, “X2”, …).
#’ @param idx Character vector of column names (or integer indices) to plot.
#’   Defaults to `c("X1","X2","X3","X4")`.  These variables will be arranged
#’   in a multi‐panel layout.
#’
#’ @return Invisibly returns `NULL`.  The function creates a multi‐panel plot of
#’   running inclusion probabilities and horizontal lines at each variable’s final
#’   (empirical) inclusion frequency.
#’
#’ @details
#’ For each variable in `idx`, the function computes
#’ \deqn{\hat p_t = \frac1t \sum_{i=1}^t v_i,}
#’ where \(v_i\in\{0,1\}\) is the indicator at iteration \(i\).  It then plots
#’ \(\hat p_t\) vs.\ iteration index, with a dashed red line at the overall
#’ \(\hat p_T\) (final mean).
#’
#’ @examples
#’ \dontrun{
#’ # Suppose `chain` is a 1000×16 matrix of 0/1 indicators named X1..X16
#’ Gibbs.RS.p.chart(chain)
#’
#’ # Or plot convergence only for predictors 5, 10, and 12
#’ Gibbs.RS.p.chart(chain, idx = c("X5","X10","X12"))
#’ }
#’
#’ @export
Gibbs.RS.p.chart = function(indicator.samples,
                            #idx=1:4,
                            idx=c("X1", "X2", "X3", "X4")){
  
  
  # mar documentation:
  # A numerical vector of the form c(bottom, left, top, right)
  # which gives the number of lines of margin to be specified on the four sides of the plot
  # The default is c(5, 4, 4, 2) + 0.1.
  
  #par(
  #  mfrow = c(ceiling(length(idx)/2), 2),
  #  mar   = c(1,1,2,1),       # almost no inner margin
  #  oma   = c(4, 4, 0, 0)     # outer margin: bottom, left, top, right
  #  )
  
  n_plots <- length(idx)
  n_rows  <- ceiling(n_plots / 2)
  
  # use more reasonable margins so xlab/ylab fit
  par(mfrow = c(n_rows, 2),
      mar   = c(4, 4, 2, 1) + 0.1,   # bottom, left, top, right
      mgp   = c(3, 1, 0))            # axis-title, axis-label, axis-line
  
  for (i in idx) {
    var_samples <- indicator.samples[, i]
    running_p <- cumsum(var_samples) / seq_along(var_samples)
    
    eps_ticks_limit=0.000005
    
    plot(
      running_p, type = 'l',
      xlab = 'Iteration',
      ylab = 'Probability',
      main = paste0("Running Inclusion Prob. Estimate for ", i),
      ylim = c(min(running_p) - eps_ticks_limit,
               max(running_p) + eps_ticks_limit),
      col = 'blue', lwd = 2
    )
    abline(h = mean(var_samples), col = 'red', lty = 2) # Final mean
  }
}


library(dplyr)

#’ Compute convergence diagnostics for MCMC inclusion indicators
#’
#’ Given a matrix of 0/1 inclusion samples from an MCMC chain (iterations × p),
#’ this function returns per‐variable metrics that quantify how stably each
#’ predictor’s inclusion probability has converged.
#’
#’ @param samples Numeric matrix with dimensions \(T \times p\), where each row
#’   is one MCMC sample of the binary inclusion indicators for all \(p\)
#’   predictors.
#’
#’ @return A tibble with one row per variable and the following columns:
#’ \describe{
#’   \item{variable}{Name of the predictor (column name of `samples`).}
#’   \item{p_full}{Final inclusion frequency \(\tfrac1T\sum_{t=1}^T v_{t,k}\).}
#’   \item{D}{Absolute difference between the first‐half and full‐chain
#’             inclusion frequencies:
#’             \(\bigl|\bar p_{\,T} - \bar p_{\,\lfloor T/2\rfloor}\bigr|\).}
#’   \item{M}{Maximum absolute deviation in the second half of the chain:
#’            \(\max_{t>\!T/2}\bigl|\hat p_t - \hat p_T\bigr|\).}
#’   \item{MCSE}{Monte Carlo standard error proxy:
#’               \(\sqrt{\hat p_T(1-\hat p_T)/T}\).}
#’ }
#’
#’ @details
#’ \itemize{
#’   \item \strong{D} flags variables whose overall inclusion probability has
#’     shifted substantially from the first half to the end of the chain.
#’   \item \strong{M} captures the worst drift late in the chain, indicating
#’     how much the running estimate can still move.
#’   \item \strong{MCSE} is a rough measure of Monte Carlo uncertainty under the
#’     assumption of i.i.d.\ Bernoulli draws.
#’ }
#’
#’ @examples
#’ \dontrun{
#’ # Suppose `chain` is a 1000×50 matrix of 0/1 inclusion indicators
#’ conv_df <- convergence_metrics(chain)
#’ # View the 10 variables with largest first‐half/final discrepancy
#’ head(conv_df, 10)
#’ }
#’
#’ @importFrom tibble tibble
#’ @importFrom dplyr arrange
#’ @export
#convergence_metrics <- function(samples) {
#  T <- nrow(samples)
#  half <- floor(T/2)

# running inclusion probabilities per variable
# seq_along(x): creates an integer sequence (1,2,..) of same length as x
# col is actual data from samples filtered by col
#  running_p <- apply(samples, 2, function(col) cumsum(col) / seq_along(col))
#  p_full   <- running_p[T, ]
#  p_half   <- running_p[half, ]

# 1. Half‐chain difference for all variables
#  D <- abs(p_full - p_half)

# 2. Max deviation in last half
# find the largest |running_p(t) - running_p(T)| for t > T/2
# for each variable k=1,..,p, 
#    for each t=T/2,..,T,
#       compute |running_p(t) - running_p(T)|
#    map k to the maximum of |running_p(t) - running_p(T)|
#  M <- sapply(seq_len(ncol(running_p)), function(k) {
#    x_k <- running_p[(half+1):T, k]
#    max(abs(x_k - p_full[k]))
#  })


# 3. Monte Carlo SE proxy
#  MCSE <- sqrt(p_full * (1 - p_full) / T)

# sort in descending order by D
#  tibble(
#    variable = colnames(samples),
#    p_full   = p_full,
#    D        = D,
#    M        = M,
#    MCSE     = MCSE
#  ) %>% arrange(desc(D))
#}


# Compute ESS via Geyer's initial positive sequence using acf()
.ess_ips <- function(x, max_lag = NULL) {
  T <- length(x)
  if (is.null(max_lag)) max_lag <- floor(T/2)
  # demean to stabilize acf
  x_cent <- x - mean(x)
  rho <- as.numeric(stats::acf(x_cent, lag.max = max_lag, plot = FALSE, type = "correlation")$acf)[-1L]
  if (length(rho) == 0) return(T)  # fall back to iid
  
  # sums of adjacent pairs: rho_{1}+rho_{2}, rho_{3}+rho_{4}, ...
  np <- length(rho) %/% 2
  pair_sum <- rho[2*(1:np)-1] + rho[2*(1:np)]
  # keep only the initial run of positive pair sums
  k <- which(pair_sum <= 0)
  if (length(k)) pair_sum <- pair_sum[seq_len(min(k)-1L)]
  
  tau <- 1 + 2*sum(pair_sum)    # integrated autocorrelation time
  if (!is.finite(tau) || tau <= 0) tau <- 1
  ESS <- T / tau
  ESS <- max(1, min(T, ESS))    # clamp to [1, T]
  ESS
}

#’ Compute convergence diagnostics for MCMC inclusion indicators (with ESS)
convergence_metrics <- function(samples, method = c("ips", "coda"), acf_max_lag = NULL) {
  method <- match.arg(method)
  T <- nrow(samples)
  half <- floor(T/2)
  
  # running inclusion probabilities per variable
  running_p <- apply(samples, 2, function(col) cumsum(col) / seq_along(col))
  p_full   <- running_p[T, ]
  p_half   <- running_p[half, ]
  
  # 1) Half‐chain difference
  D <- abs(p_full - p_half)
  
  # 2) Max deviation in last half
  M <- sapply(seq_len(ncol(running_p)), function(k) {
    x_k <- running_p[(half+1):T, k]
    max(abs(x_k - p_full[k]))
  })
  
  # 3) ESS per variable
  if (method == "coda" && requireNamespace("coda", quietly = TRUE)) {
    ESS <- as.numeric(coda::effectiveSize(samples))
  } else {
    ESS <- apply(samples, 2, .ess_ips, max_lag = acf_max_lag)
  }
  IACT <- T / ESS
  
  # MCSEs
  MCSE_iid  <- sqrt(p_full * (1 - p_full) / T)
  MCSE_ESS  <- sqrt(p_full * (1 - p_full) / ESS)
  
  tibble::tibble(
    #variable = colnames(samples) %||% paste0("V", seq_len(ncol(samples))),
    variable = colnames(samples),
    p_full   = p_full,
    D        = D,
    M        = M,
    ESS      = ESS,
    IACT     = IACT,
    MCSE_iid = MCSE_iid,
    MCSE     = MCSE_ESS   # adjusted by ESS
  ) |> dplyr::arrange(dplyr::desc(D))
}




#’ Plot p-charts for the slowest-converging variables,
#' according to some metric.
#’
#’ @param samples   MCMC indicator matrix (T × p)
#’ @param metric    One of "D", "M", or "MCSE". See convergence_metrics() doco.
#’ @param top_n     How many of the slowest variables to plot (default = 5)
#’ @export
p_chart_slowest_converging_vars <- function(samples,
                                            metric,
                                            top_n = 5) {
  # validate correctness of metric
  metric <- match.arg(metric, c("D", "M", "MCSE"))
  
  # compute convergence metrics
  cm <- convergence_metrics(samples)
  
  # pick top_n variables by descending metric
  slowest_vars <- cm %>%
    slice_max(order_by = .data[[metric]], n = top_n) %>%
    pull(variable)
  
  # plot p-charts for those
  Gibbs.RS.p.chart(samples, idx = slowest_vars)
  
  invisible(slowest_vars)
}

# ⚠️ WARNING THIS MIGHT BE OUTDATED SINCE find_BIC_minimzer() is alrd updated
# probably might delete this function since we store synthhetic datasets
# not as separate files but on a big list
format_summary_table = function(path_base=here("output", "simulation_factorial_analysis_outputs"),
                                meta=tibble(
                                  combination = 1:15,
                                  correlation = rep(c("0", "+0.5", "+0.9", "-0.5", "-0.9"), 3),
                                  s_n = rep(c("high", "moderate", "weak"), each = 5)
                                )){
  
  # Ground truth variables
  ground_truth_vars <- paste0("X", 1:4)
  
  # Initialize empty list to collect rows
  all_results <- vector("list", length = nrow(meta))
  
  for (i in 1:nrow(meta)) {
    combo_id <- paste0("c", i)
    
    # Read relative frequencies
    rel_freq_obj <- readRDS(file.path(path_base, paste0("gibbs.RS.", combo_id, ".rel.freqs.rds")))
    rel_freqs <- rel_freq_obj$rel.freqs.non.basin
    
    # Read empirical mode (most common model in MCMC)
    empirical_mode <- readRDS(file.path(path_base, paste0("gibbs.RS.", combo_id, ".best.Xs.rds")))
    
    # Read BIC-selected model
    bic_model_obj <- readRDS(file.path(path_base, paste0("gibbs.RS.", combo_id, ".best.Xs.BIC.rds")))
    bic_model <- names(which(bic_model_obj$best.model == 1))
    min_bic <- bic_model_obj$min.BIC
    
    # Logical check: were all GT variables recovered?
    recovered_all_gt_emp <- all(ground_truth_vars %in% empirical_mode)
    recovered_all_gt_bic <- all(ground_truth_vars %in% bic_model)
    
    # Combine into one row (as tibble)
    row <- tibble(
      Combination = i,
      Correlation = meta$correlation[i],
      Signal_to_Noise = meta$s_n[i],
      Empirical_Mode = paste(empirical_mode, collapse = ", "),
      BIC_Model = paste(bic_model, collapse = ", "),
      Min_BIC = min_bic,
      Recovered_GT_Empirical = recovered_all_gt_emp,
      Recovered_GT_BIC = recovered_all_gt_bic
    ) %>%
      bind_cols(as_tibble(t(rel_freqs)))  # Add relative freq columns
    
    all_results[[i]] <- row
  }
  
  # Final summary table
  return(bind_rows(all_results))
}

# install.packages(c("dplyr","knitr","kableExtra"))
# a little helper to turn one SNR‐subset into a LaTeX table
format_snr_table <- function(df, snr_level) {
  sub <- df %>%
    filter(Signal_to_Noise == snr_level) %>%
    dplyr::select( # note select() can be masked by other libraries.
      `Correlation` = Correlation,
      `v1 freq` = X1,
      `v2 freq` = X2,
      `v3 freq` = X3,
      `v4 freq` = X4,
      `Empirical mode`    = Empirical_Mode,
      `BIC model`         = BIC_Model,
      `Min. BIC`          = Min_BIC
    )
  
  cap <- paste0("Summary statistics for \\textbf{", snr_level,
                "} signal‐to‐noise across correlation levels")
  
  sub %>%
    kable(
      format       = "latex",
      booktabs     = TRUE,
      caption      = cap,
      label        = paste0("tab:", snr_level),
      align        = c("l", rep("c", ncol(sub)-1)),
      escape       = FALSE
    ) %>%
    kable_styling(latex_options = c("hold_position","scale_down")) %>%
    column_spec(1, bold = TRUE)
}

# generate result tables for latex report

format_sim_study_tables = function(summary_table_filtered){
  # now loop over each SNR level and cat() the tables
  snrs <- unique(summary_table_filtered$Signal_to_Noise)
  for (snr in snrs) {
    cat("\n\n", format_snr_table(summary_table_filtered, snr), "\n\n")
  }
}



library(dplyr)
library(knitr)
library(kableExtra)

make_summary_from_mcmc_samples <- function(mcmc.samples,
                                           ground_truth_vars = paste0("X", 1:4),
                                           data_list, ...){
  # rhos, SNRs must match your mcmc.samples keys
  
  
  rhos <- c('0', '0.5', '0.9', '-0.5', '-0.9')
  snrs <- c('strong', 'moderate', 'weak')
  
  out <- list()
  row_idx <- 1
  
  
  for(rho in rhos) for(snr in snrs){
    
    samples <- mcmc.samples[[rho]][[snr]]
    
    if (is.null(samples) || identical(samples, "empty")){
      next 
    }
    
    # data_list should be a list of lists
    # -- each element is a synthetic X and y
    # must be a data.frame
    data <- data.frame(data_list[[rho]][[snr]]$X.syn,
                       data_list[[rho]][[snr]]$y.syn)
    
    # 1. Relative inclusion freq (per column/var)
    rel_freq <- colMeans(samples)
    names(rel_freq) <- paste0("X", seq_along(rel_freq))
    
    # 2. Empirical mode (most common model)
    empirical_mode_idx <- which.max(duplicated(as.data.frame(samples)))
    empirical_mode <- names(which(samples[empirical_mode_idx, ] == 1))
    
    # 3. BIC minimizing model
    # note now find_BIC_minimizer() returns the top 5 models with smallest BIC
    
    bic_out <- find_BIC_minimizer(data, samples)
    bic_model <- bic_out$models[[1]]$selected.variables
    min_bic <- bic_out$models[[1]]$BIC
    
    # 4. GT recovery
    recovered_all_emp <- all(ground_truth_vars %in% empirical_mode)
    recovered_all_bic <- all(ground_truth_vars %in% bic_model)
    
    # Store result
    out[[row_idx]] <- tibble(
      Correlation = rho,
      Signal_to_Noise = snr,
      Empirical_Mode = paste(empirical_mode, collapse = ", "),
      BIC_Model = paste(bic_model, collapse = ", "),
      Min_BIC = min_bic,
      Recovered_GT_Empirical = recovered_all_emp,
      Recovered_GT_BIC = recovered_all_bic
    ) %>% bind_cols(as_tibble(t(rel_freq)))
    row_idx <- row_idx + 1
  }
  bind_rows(out)
}


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(coda)    # effectiveSize
  library(zoo)     # rollmean for smoothing
})

#install.packages("zoo")

# 1) Inclusion probs + ESS-based CIs -----------------------------------------
inclusion_summary <- function(samples, top_k = 10, gt_idx = NULL, ensure_gt = TRUE) {
  stopifnot(is.matrix(samples))
  #pnames <- colnames(samples) %||% paste0("X", seq_len(ncol(samples)))
  pnames <- colnames(samples) 
  
  # per-variable inclusion prob & ESS
  p_hat <- colMeans(samples)
  ess   <- apply(samples, 2, function(x) as.numeric(effectiveSize(as.mcmc(x))))
  ess[!is.finite(ess)] <- nrow(samples)  # fallback if ESS not finite
  
  se    <- sqrt(p_hat * (1 - p_hat) / ess)
  df <- tibble(
    variable = pnames,
    p_hat    = p_hat,
    ess      = ess,
    se       = se,
    lower    = pmax(0, p_hat - 1.96 * se),
    upper    = pmin(1, p_hat + 1.96 * se),
    is_gt    = if (is.null(gt_idx)) FALSE else variable %in% pnames[gt_idx]
  ) %>% arrange(desc(p_hat))
  
  # take top_k, but make sure GTs are shown if ensure_gt=TRUE
  if (!is.null(top_k)) {
    keep <- head(df$variable, top_k)
    if (ensure_gt && any(df$is_gt)) keep <- union(keep, df$variable[df$is_gt])
    df <- df %>% filter(variable %in% keep) %>% arrange(desc(p_hat))
  }
  df
}


inclusion_summary_asy_CI <- function(samples, top_k = 10, gt_idx = NULL,
                              ensure_gt = TRUE, alpha = 0.05, epsilon = 0.05) {
  stopifnot(is.matrix(samples))
  pnames <- colnames(samples)
  
  # MCMC estimates
  p_hat <- colMeans(samples)
  ess   <- apply(samples, 2, function(x) as.numeric(effectiveSize(as.mcmc(x))))
  ess[!is.finite(ess)] <- nrow(samples)
  
  # Conservative variance bound (B^2)
  var_hat <- p_hat * (1 - p_hat)
  B <- sqrt(var_hat / ess) * sqrt(nrow(samples))  # ~ n^(1/2) Var(e_n)^(1/2)
  
  z <- qnorm(1 - alpha / 2)
  
  # Asymptotic conservative CI
  half_width <- (1 + epsilon) * z * sqrt(var_hat / ess)
  lower <- pmax(0, p_hat - half_width)
  upper <- pmin(1, p_hat + half_width)
  
  df <- tibble(
    variable = pnames,
    p_hat    = p_hat,
    ess      = ess,
    lower    = lower,
    upper    = upper,
    is_gt    = if (is.null(gt_idx)) FALSE else variable %in% pnames[gt_idx]
  ) %>% arrange(desc(p_hat))
  
  if (!is.null(top_k)) {
    keep <- head(df$variable, top_k)
    if (ensure_gt && any(df$is_gt)) keep <- union(keep, df$variable[df$is_gt])
    df <- df %>% filter(variable %in% keep) %>% arrange(desc(p_hat))
  }
  
  df
}




plot_inclusion_custom_fill <- function(df,
                           title = "Inclusion probabilities (w/ ESS 95% CI)",
                           col_fill=NaN) {
  
  if (!is.na(col_fill)){
    df = df %>% mutate(col_fill=factor(.data[[col_fill]],
                                       levels=unique(df[,col_fill])))
    
    base = ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat,
                          color = col_fill, fill=col_fill)) 
  } else{
    base = ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat)) 
  }

  base = base +
    geom_col(width = 0.7, alpha = 0.2, color = NA) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, color = "grey20") +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
    scale_color_discrete() + 
    #scale_fill_discrete() + 
    #scale_color_manual(values = c(`TRUE` = "#dc2626", `FALSE` = "#374151"), guide = "none") +
    #scale_fill_manual(values  = c(`TRUE` = "#dc2626", `FALSE` = "#9ca3af"), guide = "none") +
    labs(x = NULL, y = "Inclusion probability", title = title,
         subtitle = "Error bars use ESS-adjusted standard errors") +
    coord_flip() +
    theme_minimal(base_size = 12)
  
  base 
}


plot_inclusion <- function(df,
                           title = "Inclusion probabilities (w/ ESS 95% CI)",
                           elbow_clusters=FALSE,
                           errorbar=TRUE) {
  
  if (elbow_clusters == TRUE){
    
    df = df %>% mutate(clusters=factor(clusters, levels=unique(df$clusters)))
    
    base = ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat,
                          color = clusters, fill=clusters)) 
    base = base +
      geom_col(width = 0.7, alpha = 0.2, color = NA) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, color = "grey20") +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
      scale_color_discrete() + 
      #scale_fill_discrete() + 
      #scale_color_manual(values = c(`TRUE` = "#dc2626", `FALSE` = "#374151"), guide = "none") +
      #scale_fill_manual(values  = c(`TRUE` = "#dc2626", `FALSE` = "#9ca3af"), guide = "none") +
      labs(x = NULL, y = "Inclusion probability", title = title,
           subtitle = "Error bars use ESS-adjusted standard errors") +
      coord_flip() +
      theme_minimal(base_size = 12)
  } else{
    base = ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat,
                          color = is_gt, fill = is_gt)) +
      geom_col(width = 0.7, alpha = 0.2, color = NA)
    
    if (errorbar){
      base = base + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25,
                           color = "grey20", alpha=0.4)
    }
      
    base = base +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
      scale_color_manual(values = c(`TRUE` = "#dc2626", `FALSE` = "#374151"), guide = "none") +
      scale_fill_manual(values  = c(`TRUE` = "#dc2626", `FALSE` = "#9ca3af"), guide = "none") +
      labs(x = NULL, y = "Inclusion probability", title = title,
           subtitle = "Error bars use ESS-adjusted standard errors") +
      coord_flip() +
      theme_minimal(base_size = 12)
  }
  
  base 
}


# fixed diverging palette: red/orange (<0.5)  →  neutral  →  teal/blue (>0.5)
prob_scales <- function() {
  list(
    scale_color_gradient2(
      low = "#d73027",   # red
      mid = "#fee08b",   # warm neutral
      high = "#1f9e89",  # teal/blue
      midpoint = 0.5,
      limits = c(0, 1),
      oob = scales::squish,
      guide = "none"
    ),
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee08b",
      high = "#1f9e89",
      midpoint = 0.5,
      limits = c(0, 1),
      oob = scales::squish,
      guide = "none"
    )
  )
}

# Make a discrete hue palette for any number of clusters
mk_cluster_palette <- function(levels, base = "Set3") {
  n <- length(levels)
  base_cols <- RColorBrewer::brewer.pal(min(max(3, n), 12), base)
  # extend smoothly if n > length(base palette)
  cols <- if (n <= length(base_cols)) base_cols[seq_len(n)]
  else grDevices::colorRampPalette(base_cols)(n)
  setNames(cols, levels)
}


plot_inclusion_v2 <- function(df,
                           title = "Inclusion probabilities (w/ ESS 95% CI)",
                           color_by = c("prob","clusters"),
                           cluster_levels = NULL,      # optional: lock levels order
                           cluster_palette = NULL) {   # optional: pass a fixed palette
  
  color_by <- match.arg(color_by)
  
  if (color_by == "clusters") {
    # fix levels so colors are consistent across plots
    if (is.null(cluster_levels)) cluster_levels <- sort(unique(df$clusters))
    df <- df %>% dplyr::mutate(clusters = factor(clusters, levels = cluster_levels))
    
    # build or use provided palette
    if (is.null(cluster_palette)) cluster_palette <- mk_cluster_palette(levels(df$clusters))
    #print(df)
    if (max(df$p_hat) == 1){
      UPPER_LIMIT = 1
    } else {
      UPPER_LIMIT = max(df$upper, na.rm = TRUE) + 0.01
    }
    
    p <- ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat,
                        color = clusters, fill = clusters)) +
      #ylim(0, max(df$upper) + 0.1) + 
      geom_col(width = 0.7, alpha = 0.20, color = NA) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, color = "grey25") +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0,UPPER_LIMIT),
                         expand = expansion(mult = c(0.02, 0.05))) +
      scale_fill_manual(values = cluster_palette, drop = FALSE) +
      scale_color_manual(values = cluster_palette, drop = FALSE, guide = "none") +
      coord_flip() +
      labs(x = NULL, y = "Inclusion probability", title = title,
           subtitle = "Error bars use ESS-adjusted standard errors") +
      theme_minimal(base_size = 12)
    
  } else {
    # (keeps the continuous prob palette option if you still want it)
    p <- ggplot(df, aes(x = reorder(variable, p_hat), y = p_hat,
                        color = p_hat, fill = p_hat)) +
      geom_col(width = 0.7, alpha = 0.20, color = NA) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, color = "grey25") +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0.02, 0.05))) +
      scale_color_gradient2(low="#d73027", mid="#fee08b", high="#1f9e89",
                            midpoint=0.5, limits=c(0,1), oob=scales::squish, guide="none") +
      scale_fill_gradient2(low="#d73027", mid="#fee08b", high="#1f9e89",
                           midpoint=0.5, limits=c(0,1), oob=scales::squish, guide="none") +
      coord_flip() +
      labs(x = NULL, y = "Inclusion probability", title = title,
           subtitle = "Error bars use ESS-adjusted standard errors") +
      theme_minimal(base_size = 12)
  }
  
  p
}





# 2) Trace of number of active variables -------------------------------------
plot_num_vars_trace <- function(samples, window = 50,
                                title = "Trace: number of variables in sampled models",
                                dimension_BIC_minim = 5) {
  # compute number of variables in each sampled model, as well as rolling mean
  nv   <- rowSums(samples)
  smth <- zoo::rollmean(nv, k = window, fill = NA, align = "right")
  df   <- tibble(iter = seq_along(nv),
                 n_vars = nv,
                 n_vars_smooth = smth)
  
  ggplot(df, aes(iter, n_vars)) +
    
    # raw trace
    geom_line(alpha = 1, color = "black") +
    
    # rolling mean with legend
    geom_line(aes(y = n_vars_smooth, color = "Rolling mean"),
              linewidth = 1.1, na.rm = TRUE) +
    
    # horizontal BIC minimizer line with legend
    geom_hline(aes(yintercept = dimension_BIC_minim, color = "BIC minimizer"), linewidth = 1) +
    
    # annotation for red line
    annotate("text", x = -10, y = dimension_BIC_minim - 0.7,
             label = dimension_BIC_minim, color = "red", hjust = 0) +
    
    labs(x = "Iteration", y = "# variables in model", title = title, color = "") +
    scale_color_manual(values = c("Rolling mean" = "#2563eb",
                                  "BIC minimizer" = "red")) +
    
    theme_minimal(base_size = 12) +
    ylim(0, max(nv) + 2) +
    
    # note legend is automatically generated whenever there is aes()
    theme(legend.position = "top",
          legend.title = element_blank()) +
    
    # remove grid lines
    #theme(panel.grid.major = element_blank(),
    #      panel.grid.minor = element_blank()) 
    
    theme(panel.grid.minor = element_blank())
}


# 3) Heatmap of model inclusion over iterations -------------------------------
# note we can tweak colors using hsv(hue, saturation, value)
plot_model_heatmap <- function(
    samples,
    gt_idx = 1:4,
    thin = 1,
    order = c("by_mean", "as_is"),
    title = "Model path heatmap",
    exclude_unimportant_vars = 0,
    put_gts_on_top = TRUE,
    transpose = FALSE,
    colors = list(
      GT_1     = "firebrick",
      NonGT_1  = "steelblue",
      EXCLUDED = "gray90"
    )
    
) {
  order <- match.arg(order)
  stopifnot(is.matrix(samples) || is.data.frame(samples))
  samples <- as.matrix(samples)
  Tn <- nrow(samples); p <- ncol(samples)
  pnames <- colnames(samples) 
  
  ## thinning
  keep_iter <- seq(1, Tn, by = max(1, thin))
  S <- samples[keep_iter, , drop = FALSE]
  
  ## order variables
  ord <- seq_len(p)
  if (order == "by_mean") {
    score <- colMeans(S)
    ord <- order(score, decreasing = TRUE)
    if (put_gts_on_top && !is.null(gt_idx)) {
      gt_names <- pnames[gt_idx]
      ord <- c(match(gt_names, pnames), setdiff(ord, match(gt_names, pnames)))
      ord <- unique(ord[!is.na(ord)])
    }
  }
  
  ## keep only top variables after pruning
  keep_p <- max(1, min(length(ord), length(ord) - exclude_unimportant_vars))
  ord <- ord[seq_len(keep_p)]
  S   <- S[, ord, drop = FALSE]
  vp  <- pnames[ord]
  
  ## long data
  long <- as.data.frame(S)
  long$iter <- keep_iter
  long <- tidyr::pivot_longer(long, -iter, names_to = "variable", values_to = "in_model")
  long$variable <- factor(long$variable, levels = rev(vp))
  long$is_gt <- ifelse(long$variable %in% pnames[gt_idx], "GT", "NonGT")
  
  ## composite fill
  long$fill_state <- ifelse(long$in_model == 1,
                            paste0(long$is_gt, "_1"),
                            "EXCLUDED")
  long$fill_state <- factor(long$fill_state,
                            levels = c("GT_1","NonGT_1","EXCLUDED"))
  
  ## palette + labels
  pal <- c("GT_1"     = colors$GT_1,
           "NonGT_1"  = colors$NonGT_1,
           "EXCLUDED" = colors$EXCLUDED)
  labs_fill <- c("GT_1" = "Ground truth included",
                 "NonGT_1" = "Non-GT included",
                 "EXCLUDED" = "Excluded")
  
  ## base plot
  if (!transpose) {
    g <- ggplot(long, aes(x = iter, y = variable, fill = fill_state)) +
      geom_raster() +
      scale_fill_manual(values = pal, labels = labs_fill, name = "State") +
      labs(x = "Iteration", y = NULL, title = title) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.margin = unit(c(0.00001, 0.1, 0.00001, 0.1), "cm"))
  } else {
    g <- ggplot(long, aes(x = variable, y = iter, fill = fill_state)) +
      geom_raster() +
      scale_fill_manual(values = pal, labels = labs_fill, name = "State") +
      labs(x = NULL, y = "Iteration", title = title) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.margin = unit(c(0.00001, 0.1, 0.00001, 0.1), "cm"))
  }
  
  g
}





