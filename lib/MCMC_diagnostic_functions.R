#' Given multiple runs of MCMC samples, this function plots the interleaved
#' trace plots of variables in a grid (e.g. top 5-6) with the highest positive or
#' negative posterior mean (*), since plotting trace plots
#' of all variables in a grid wouldn't be feasible. 
#'
#' (*) NOTE: the top variables are selected based on posterior means computed
#' from the FIRST chain only
#'
#'
#' @param ... a variable number of samples generated from an MCMC algorithm,
#' each being a [iter x p] matrix, where p is the dimension of beta coeffs
#'
#' @param optional 'positive' or 'negative' -- whether the user wants to look at 
#' trace plots of coefficients with either positive largest or negative largest
#' posterior means. Defaults to 'positive'.
#'
#' @param CUTOFF an optional vertical line to be put on the respective 
#' x-axis value. Useful to detect the need for extending the burn-in period.
#' Defaults to NaN.
#'
#' @return Produces trace plots for the top 6 coefficients (based on the 1st chain)
#' 
#' @examples
#' traceplots.high.posterior.means(chain1, chain2, sign="negative", CUTOFF=1000)


traceplots.high.posterior.means = function(...,
                                           sign='positive',
                                           CUTOFF=NaN,
                                           save_png = FALSE,
                                           file = "trace_plots.png"){
  
  # 1. Consolidating the variable number of chains given.
  chains = list(...)
  
  if (length(chains) == 0) {
    stop("Please provide at least a MCMC chain in the form of a matrix.")
  }
  
  ## Posterior mean selection for variable ranking is based on first chain
  reference_chain = chains[[1]]
  
  # 2. Get top 6 variables with largest posterior means in the chosen sign/
  # direction
  
  if (sign == 'negative'){
    top.feature.coeffs = names(sort(colMeans(reference_chain), decreasing=FALSE))[1:6]
  } else{
    top.feature.coeffs = names(sort(colMeans(reference_chain), decreasing=TRUE))[1:6]
  }
  
  # 3. Begin plotting
  
  if (save_png) png(file, width = 1500, height = 800)
  
  par(mfrow = c(2, 3))
  
  
  for (top.feature in top.feature.coeffs){
    
    ## 3.1 Plot the first chain
    
    # posterior mean of reference chain
    posterior.mean.reference = round(mean(reference_chain[, top.feature]), 4)
    
    plot(1:nrow(reference_chain),
         reference_chain[, top.feature],
         type = "l", col = 1,
         main = paste0("Trace plot for\n", top.feature),
         xlab = paste0("Iteration\n", 
                       "Posterior mean (reference chain): ",
                       posterior.mean.reference),
         ylab = "Value")
    
    
    ## 3.2 Overlay remaining trace plots, if more than 1 chains are supplied
    
    if (length(chains) > 1) {
      for (k in 2:length(chains)) {
        
        lines(1:nrow(chains[[k]]),
              chains[[k]][, top.feature],
              col = k)
      }
    }
    
    ## 3.3 Burn-in reference line. 
    if (!is.nan(CUTOFF)){
      abline(v=CUTOFF, col='red', lty=2)
    }
    
    
    # Posterior mean label (from first chain)
    #mtext(paste("Posterior mean =", round(mean(reference_chain[, top.feature]), 4)),
    #      side = 1, line = 3, cex = 0.8)
  }
  
  if (save_png) dev.off()
  
}



plot_acf.values_intervals = function(chain, lags = 50,
                              save_png = FALSE,
                              file = "acf_intervals.png") {
  
  # 1. Initialize data structures
  p = ncol(chain)
  variable_names = colnames(chain)
  
  # 2 vectors to store the endpoints of the acfs
  min_acf = numeric(p)
  max_acf = numeric(p)
  
  
  # 2. Obtain the min and max of the acf values per variable
  
  for (j in 1:p) {
    
    acf_vals = acf(chain[, j], plot = FALSE, lag.max = lags)$acf
    acf_lags = acf_vals[-1]  # exclude lag 0
    min_acf[j] = min(acf_lags)
    max_acf[j] = max(acf_lags)
    
  }
  
  if (save_png) png(file, width = 1500, height = 800)
  
  # 3. set up a skeleton of the plot 
  
  plot(1:p, rep(NA, p), ylim = c(min(min_acf), max(max_acf)),
       xaxt = "n", xlab='', ylab = "ACF (lag > 0)",
       main = "ACF Interval Ranges per Variable")
  
  # set up line segments to connect (i,min_acf) to (i,max_acf)
  # i.e. min_acf--------max_acf for all variables i
  segments(1:p, min_acf, 1:p, max_acf, col = "steelblue", lwd = 2)
  
  # emphasize max, min of acf values by labelling them
  points(1:p, max_acf, pch = 16, col = "darkred")  # top ACF per variable
  points(1:p, min_acf, pch = 16, col = "darkgreen")  # min ACF per variable
  
  
  axis(1, at = 1:p, labels = variable_names, las = 2, cex.axis = 0.6)
  abline(h = 0, col = "gray60", lty = 2)
  
  if (save_png) dev.off()
}


library(coda)
gelman.rubin.diagnostics = function(...){
  
  # 1. Consolidating the variable number of chains given.
  chains = list(...)
  
  if (length(chains) == 0) {
    stop("Please provide at least a MCMC chain in the form of a matrix.")
  }
  
  # 2. in the case where the chains have different length, set all
  # to be the minimum one
  
  agreed.chain.length = min(sapply(chains,nrow))
  
  # 2.1 match all chains to have the same (min) length
  chains.trimmed = lapply(chains,
                          function(chain)chain[1:agreed.chain.length, , drop=FALSE])
  
  # 3. need to convert data matrix to a mcmc.list format
  mcmc.chains = mcmc.list(lapply(chains.trimmed, mcmc))
  
  # 4. do the test -- compute the GR test statistic
  return(gelman.diag(mcmc.chains, autoburnin = FALSE))
  
}





