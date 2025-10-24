#' Obtain samples from the posterior of beta and latent variable
#' for Probit logistic regression with ridge prior with a fixed hyperparameter
#' (lambda)
#'
#'
#' @param X = the n x p design matrix, MUST include intercept column
#' @param y = the response vector (a nx1 matrix)
#' @param iter = number of iterations in the Gibbs sampler
#' @param burnin = number of burn-in samples
#' @param beta0 = initial values for the beta vectors. MUST be a px1 matrix.
#' @param lambda = the ridge hyperparameter (fixed) 
#'
#' @returns a matrix of size (iter-burnin) x (p), where p is the 
#' dimension of beta 

Probit.Gibbs.Sampler = function(X, y, iter, burnin, beta0, lambda){
  
  start.time.function = Sys.time()
  
  # 0. precalculate quantities for parameters of conditional dist
  p = ncol(X); n=nrow(X)
  XTX = crossprod(X)
  ridge.XTX.inv = solve( (XTX + lambda*diag(p)) )
  ridge.hat.matrix = ridge.XTX.inv %*% t(X) 
  samples = matrix(0, iter, p); colnames(samples) = colnames(X)
  beta = as.matrix(beta0)
  # no need for z=as.matrix(numeric(n), nrow=n, ncol=1),
  # as R does allow multiplying a matrix by a vector
  z=numeric(n)
  
  # for truncated normal samplings
  idx.pos = which(y == 1)
  idx.neg = which(y == 0)
  
  # due to numerical precision issues, we have 
  # to force the covariance matrix for cond. posterior of beta
  # to be symmetric.
  ridge.XTX.inv.corrected = 0.5 * (ridge.XTX.inv + t(ridge.XTX.inv))
  
  # 1. begin the sampling process
  
  # track percentage of completion in increments of 10% of iter
  PROP.COMPLETE = 0
  
  for (i in 1:iter){
    
    if (i >= PROP.COMPLETE * iter){
      cat(100*PROP.COMPLETE,"% completed",
          " i.e. ", PROP.COMPLETE * iter, "iterations",
          "\n")
      PROP.COMPLETE = PROP.COMPLETE + 0.1
    } 
    
    
    # 1.1 sample z_1 ,.., z_n using its conditional distribution z_i^(i)|beta^(i-1), y
    
    # instead of looping through each z_i, use vectorized operations to 
    # sample in bulk for each binary category (y=1 and 0) and assign directly
    
    # for(j in 1:n){
    #  mu.j = X[j,] %*% beta
    #  
    #  if( y[j]==1 ){
    #    z[j] = rtruncnorm(1, a=0, b=Inf, mean=mu.j, sd=1)
    #  } else {
    #    z[j] = rtruncnorm(1, a=-Inf, b=0, mean=mu.j, sd=1)
    #  }
    #}
    
    # drop=FALSE is to keep it as a matrix/vector in case idx_pos or idx_neg
    # is only 1 element
    
    #cat("dim of X", dim(X), '\n')
    #cat("dim of beta", dim(beta) , '\n')
    
    mu = X %*% beta
    z[idx.pos] = rtruncnorm(length(idx.pos),
                            a=0, b=Inf, mean=mu[idx.pos, , drop=FALSE], sd=1)
    
    z[idx.neg] = rtruncnorm(length(idx.neg),
                            a=-Inf, b=0, mean=mu[idx.neg, , drop=FALSE], sd=1)
    
    # 1.2 sample beta from its conditional distribution beta^(i)|z_i^(i)
    beta = rmvnorm(1,
                   mean=as.vector(ridge.hat.matrix %*% z),
                   sigma=ridge.XTX.inv.corrected )
    
    # result is a row vector, need to reshape it into a col matrix
    # for future matrix operations.
    beta = t(beta)
    
    
    # 1.3 append the sampled beta
    samples[i,] = beta
    
    
  }
  
  end.time.function = Sys.time()
  cat("\n")
  time.diff = end.time.function - start.time.function
  cat("Overall time taken (minutes): ", as.numeric(time.diff, units = "mins"))
  
  # remove first few samples as burnin
  return(samples[-c(1:burnin), ])
  
}