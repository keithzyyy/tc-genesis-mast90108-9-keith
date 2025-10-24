set.seed(123)
# one hot encoding with intercept
source(here("lib", "preprocessing_v3.R"))

# full rank encoding with NO intercept
source(here("lib", "preprocessing_v2.R"))



#' Stratified Sampling by Basin with Minimum Coverage for a Rare Class
#'
#' This function performs stratified sampling from a dataset with a 'basin' column,
#' ensuring proportional representation from each basin *except* for a specified
#' rare basin (e.g., "SA"), which is included in full.
#'
#' The remaining sampling quota (N.SIM - n_min_basin) is distributed proportionally
#' across the other basins based on their original frequencies.
#'
#' @param dta A data frame containing a categorical 'basin' variable. It is assumed
#'        that missing values may have been encoded as a literal "NA" string, which
#'        is corrected internally using `fix_NA_value_and_refactor_basin()`.
#' @param N.SIM Integer. Total number of samples to return (including the full rare basin).
#'        Default is 2500.
#' @param min.basin String. Name of the rare basin to include in full. Default is "SA".
#'
#' @return A stratified sampled data frame with `N.SIM` rows. The rare basin (`min.basin`)
#'         is included completely, and the remaining samples are proportionally drawn from
#'         other basins.
#'
#' @details
#' The steps are:
#' 1. Correct improperly coded missing values in the `basin` variable.
#' 2. Identify the number of rows in the rare basin (`min.basin`) and reserve them all.
#' 3. Distribute the remaining quota proportionally across other basins.
#' 4. Fix rounding discrepancies (due to `round`) so that the final total is exactly `N.SIM`.
#' 5. Sample each basin accordingly, using a fixed seed for reproducibility.
stratified_sample_by_basin <- function(dta,
                                       N.SIM = 2500,
                                       min.basin = "SA") {
  
  # Step 0: fix the incorrectly interpreted NA as an actual category, "NA"
  dta = fix_NA_value_and_refactor_basin(dta)
  
  # Step 1: Get counts for each basin
  basin_counts <- table(dta$basin)
  basin_levels <- names(basin_counts)
  
  # Step 2: Fix the small basin first
  n_min_basin <- basin_counts[min.basin]
  idx_min_basin <- which(dta$basin == min.basin)
  
  # Step 3: Remaining N.SIM to distribute
  remaining_N <- N.SIM - n_min_basin
  other_basins <- setdiff(basin_levels, min.basin)
  
  # Step 4: Calculate proportional allocation
  counts_other <- basin_counts[other_basins]
  prop_other <- counts_other / sum(counts_other)
  allocation <- round(prop_other * remaining_N)
  
  # Fix rounding errors to make total sum correct
  diff <- remaining_N - sum(allocation)
  if (diff != 0) {
    # Add/subtract diff to the largest/smallest group to fix
    biggest <- which.max(allocation)
    allocation[biggest] <- allocation[biggest] + diff
  }
  
  # Step 5: Perform sampling per basin
  set.seed(123)
  idx_sampled <- idx_min_basin
  
  for (b in other_basins) {
    rows_in_b <- which(dta$basin == b)
    size_b <- min(allocation[b], length(rows_in_b))  # just in case
    idx_sampled <- c(idx_sampled, sample(rows_in_b, size = size_b))
  }
  
  # Return stratified sampled raw data
  dta_sampled <- dta[sort(idx_sampled), ]
  return(dta_sampled)
}


#' Given a sample of the TC genesis data, preprocess it into
#' a design matrix by encoding categorical variables by a method
#' of the user's choice (via the 'encoding' parameter)
#'
preproc_stratified_data = function(dta_sampled, encoding){
  
  # one hot encoding with intercept
  if (encoding == "OH"){
    data.processed = preprocess_v3(dta_sampled, fix.NA=FALSE)
  }
  
  # full rank encoding without intercept
  if (encoding == "FR"){
    data.processed = preprocess_v2(dta_sampled, fix.NA=FALSE)
  }
  
  
  X.sim <- data.processed$predictors
  
  return(X.sim)
}


#' Automates the proper stratified sampling process from
#' TC_sim_regenerate_data_case_1.Rmd.
#'
sample.rows.simulation.v2 = function(data, min.basin="SA", N.SIM=1000){
  
  set.seed(123)
  # sample rows 
  sample.class.1 <- stratified_sample_by_basin(data[data['TC_genesis'] == 1,],
                                              N.SIM = N.SIM,
                                              min.basin = min.basin)
  
  sample.class.0 <- stratified_sample_by_basin(data[data['TC_genesis'] == 0,],
                                              N.SIM = N.SIM,
                                              min.basin = min.basin)
  
  
  # combine rows together
  data.simulated.samples <- rbind(sample.class.0, sample.class.1)
  
  # turn into a design matrix
  X.sim.balanced.class <- preproc_stratified_data(data.simulated.samples,
                                                 encoding="FR") 
  
  return(X.sim.balanced.class)
}




set.seed(123)

#' computes AIC of the logistic regression model comprised
#' of predictors indexed by indicator variable v
#'
#' @param data = the data frame comprised of cbind(X,y) where
#' X(nxp) is the data matrix and y(nx1) is the response. 
#' X should already include an intercept. 
#' 
#' @param v = the current state in the markov chain.
#' 
#' @param method = either use glm() or manually compute the AIC.
#' 
#' @returns AIC of the model indexed by v 
compute.AIC = function(data,v){
  
  # there may be some components of v that are 0,
  # meaning we should exclude them from the model.
  # note data is (n x (p+1)) due to response being in 
  # the last column
  X.sub = data[, which(v == 1), drop=FALSE]
  y.sub = data[, ncol(data)]
  data.model = data.frame(X.sub, y = y.sub)
  
  # Check for degenerate case
  if (ncol(X.sub) == 0 || length(unique(y.sub)) == 1) return(Inf)
  
  # note X already includes the intercept
  fit = tryCatch({
    mod = glm(y ~ 0 + ., data = data.model, family = binomial)
  }, error = function(e) return(NULL))
  
  if (is.null(fit) || is.null(fit$coefficients) || any(is.na(fit$coefficients))) {
    return(Inf)
  }
  
  # Compute probabilities
  eta = as.vector(as.matrix(X.sub) %*% fit$coefficients)
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
simulate.conditional.vk = function(k, data, v, penalty){
  
  # Want: P(v) / [P(v0) + P(v1)]
  # v0 := v with v_k = 0
  # v1 := v with v_k = 1
  v0 <- v; v0[k] <- 0
  v1 <- v; v1[k] <- 1
  
  # calculate kernel of P(v)
  #prob.joint = calculate.kernel.prob.of.model(data, v, penalty)
  
  # calculate kernel of P(v) but v[k] = 0
  #v.copy[k] = 0
  #prob.joint.0 = calculate.kernel.prob.of.model(data, v.copy, penalty)
  
  # calculate kernel of P(v) but v[k] = 1
  #v.copy[k] = 1
  #prob.joint.1 = calculate.kernel.prob.of.model(data, v.copy, penalty)
  
  # obtain the success probability
  #p.succ = prob.joint / (prob.joint.0 + prob.joint.1)
  
  if (penalty == 'AIC'){
    log.p0 = -compute.AIC(data, v0)
    log.p1 = -compute.AIC(data, v1)
  }
  
  if (penalty == 'BIC'){
    log.p0 = -compute.BIC(data, v0)
    log.p1 = -compute.BIC(data, v1)
  }
  
  
  # Use log-sum-exp trick to avoid underflow
  max.log.p = max(log.p0, log.p1)
  denom = log(exp(log.p0 - max.log.p) + exp(log.p1 - max.log.p)) + max.log.p
  
  log.p.succ = log.p0 - denom
  p.succ = exp(log.p.succ)
  
  # cat("p.succ=", p.succ, "\n")
  
  if (p.succ == 0){
    return(0)
  }
  
  
  # finally simulate from associated bernoulli dist
  return(rbinom(1,1,p.succ))
  
}


#' Perform a MCMC random search procedure to find the model
#' that minimizes some model selection criterion.
#'
#'
#' @param X = the n x p design matrix, MUST include intercept column
#' @param y = the response vector (a nx1 matrix)
#' @param v0 = initial values for p-dimensional indicator vectors
#'
#' @returns a matrix of size (iter-burnin) x (p), where p is the 
#' dimension of beta 

Gibbs.random.search = function(v0, data,
                               burnin=200,
                               iter=800,
                               penalty='AIC'){
  set.seed(123)
  
  start.time.function = Sys.time()
  
  # 0. defining essential quantities
  p = ncol(data) - 1
  samples = matrix(0, iter, p)
  samples[1,] = v0
  
  
  # 1. begin the gibbs sampler iteration!
  for (j in 2:iter){ # j = MCMC iteration idx
    
    # 1.1 obtain the most recent state of the markov chain
    # v = (v^(j-1)_1, v^(j-1)_2, ... , v^(j-1)_p)
    v.prev = samples[j-1,]
    
    # initialize the next state, which will be updated component wise
    v = v.prev
    
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
    
    
    # for diagnose purpose
    #if (j == 10){
    #  return(NaN)
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








