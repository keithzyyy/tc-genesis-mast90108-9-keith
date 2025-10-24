#!/usr/bin/env Rscript

## set metadata 
SEED=123
ITER   <- 2000
BURNIN <- 200
GLOBAL_START_TIME <- Sys.time()

set.seed(SEED) # NOTE: NO further set.seed() needs to be invoked.


## Load all external libraries 
if (!requireNamespace("here", quietly=TRUE)) install.packages("here")
library(here) # to override here(), simply go to File -> New Project
here()
library(mvtnorm)
library(truncnorm)
library(Matrix)
library(MASS)
library(mvnTest)
library(dplyr)
library(knitr)
library(kableExtra)
if (!require(future)) install.packages("future", dependencies=TRUE)
library(future)
plan(multisession, workers = 3) # set 3 cores for now


## Load all local libraries

# one hot encoding with intercept
source(here("lib", "preprocessing_v3.R"))
# full rank encoding with NO intercept
source(here("lib", "preprocessing_v2.R"))
source(here("lib", "probit_gibbs.R"))
#source(here("lib", "MCMC_diagnostic_functions.R"))
source(here("lib", "simulation_routines.R"))
source(here("lib", "factorial_analysis_routines.R"))
source(here("lib", "gibbs_sampler_random_search.R"))


## Base directory to save MCMC results to
BASE_DIR = "simulation_factorial_analysis_n_2000_p_160_outputs"


## Read synthetic datasets for the 2 studies

results_n_2000_p_160 = readRDS(
  here("output",
       "simulation_factorial_analysis_outputs",
       "results_n_2000_p_160.rds"))

X.no.corr = results_n_2000_p_160[['0']][['moderate']]$X.syn
y.no.corr = results_n_2000_p_160[['0']][['moderate']]$y.syn

X.high.corr = results_n_2000_p_160[['0.9']][['moderate']]$X.syn
y.high.corr = results_n_2000_p_160[['0.9']][['moderate']]$y.syn



## HELPER FUNCTION FOR RS RUNNING


run_three_RS = function(X, y, iter=2200, burnin=200, SEED){
  
  # although we wanna ensure v02 and v03 are diferent
  # yet reproducible
  # we do not need to set.seed(SEED) because we 
  # have already done that at the beginning of this script. O/w
  # v02 and v03 will yield the same initial values 
  
  ## init values
  p = ncol(X)
  v01 <- rep(1, p)
  v02 <- sample(0:1, p, replace = TRUE)
  v03 <- sample(0:1, p, replace = TRUE)
  
  cat("Running 3 MCMC Random Searches now. \n")
  
  
  # run multiple chains with different init values
  # and track time of each chain. Add exception handling to
  # not kill the entire funcion execution.
  
  # 1st init val
  f1 <- future({
    tryCatch({
      cat("Running Chain with 1st Initial Value.. \n")
      t0 <- Sys.time()
      chain <- suppressWarnings(
        Gibbs.random.search(v01, data.frame(X, y),
                            burnin=burnin, iter=iter, penalty='BIC')
      )
      t1 <- Sys.time()
      runtime = as.numeric(difftime(t1, t0, units="mins"))
      cat("Runtime for Chain with 1st Initial Value:", runtime, "\n")
      list(chain=chain,
           runtime=runtime)
    }, error = function(e) list(error = e$message))
  }, seed = SEED + 1)
  
  # 2nd init val
  f2 <- future({
    tryCatch({
      cat("Running Chain with 2nd Initial Value.. \n")
      t0 <- Sys.time()
      chain <- suppressWarnings(
        Gibbs.random.search(v02, data.frame(X, y),
                            burnin=burnin, iter=iter, penalty='BIC')
      )
      t1 <- Sys.time()
      runtime = as.numeric(difftime(t1, t0, units="mins"))
      cat("Runtime for Chain with 2nd Initial Value:", runtime, "\n")
      list(chain=chain,
           runtime=runtime)
    }, error = function(e) list(error = e$message))
  }, seed = SEED + 2)
  
  # 3rd init val
  f3 <- future({
    tryCatch({
      cat("Running Chain with 3rd Initial Value.. \n")
      t0 <- Sys.time()
      chain <- suppressWarnings(
        Gibbs.random.search(v03, data.frame(X, y),
                            burnin=burnin, iter=iter, penalty='BIC')
      )
      t1 <- Sys.time()
      runtime = as.numeric(difftime(t1, t0, units="mins"))
      cat("Runtime for Chain with 3rd Initial Value:", runtime, "\n")
      list(chain=chain,
           runtime=runtime)
    }, error = function(e) list(error = e$message))
  }, seed = SEED + 3)
  
  # collect results
  res1 <- value(f1)
  res2 <- value(f2)
  res3 <- value(f3)
  
  # return both chains and runtimes
  return(
    list(
      init1 = res1$chain, time1 = res1$runtime, initval1 = v01,
      init2 = res2$chain, time2 = res2$runtime, initval2 = v02,
      init3 = res3$chain, time3 = res3$runtime, initval3 = v03
    )
  )
  
  
}

## check directory

cat(sprintf("[%s] Verifying existence of output directory...\n", Sys.time()))

outdir <- here("output", BASE_DIR)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat(sprintf("[%s] Output Directory Successfully Verified. \n", Sys.time()))


## commence the 2 studies.

studies <- list(
  no_corr  = list(X = X.no.corr,  y = y.no.corr),
  high_corr= list(X = X.high.corr,y = y.high.corr)
)


for(name in names(studies)){
  tryCatch({
    
    cat(sprintf("[%s] Begin running RS for %s\n", Sys.time(), name))
    res.RS = run_three_RS(studies[[name]]$X, studies[[name]]$y, 
                          iter=ITER, burnin=BURNIN, SEED=SEED)
    
    filename <- sprintf("gibbs_3chains_results_%s.rds", name)
    saveRDS(res.RS, file.path(outdir, filename))
    cat(sprintf("[%s] Saved %s\n", Sys.time(), filename))
    
  }, error = function(e){
    cat(sprintf("[%s] ERROR on %s: %s\n", Sys.time(), name, e$message))
  })
}




## Save metadata for reproducibility
GLOBAL_END_TIME = Sys.time()
GLOBAL_TIME_TAKEN = as.numeric(difftime(GLOBAL_END_TIME,
                                        GLOBAL_START_TIME, units="mins"))

metadata <- list(
  SEED             = SEED,
  ITER             = ITER,
  BURNIN           = BURNIN,
  timestamp_start  = as.character(GLOBAL_START_TIME),
  timestamp_end    = as.character(GLOBAL_END_TIME),
  global_time_mins = GLOBAL_TIME_TAKEN,
  session          = sessionInfo()
)


cat(sprintf("[%s] Job metadata:\n", Sys.time()))
print(metadata)


cat(sprintf("[%s] Saving the job metadata.. \n", Sys.time()))
saveRDS(metadata,
        file = file.path(outdir, "metadata.rds"))










