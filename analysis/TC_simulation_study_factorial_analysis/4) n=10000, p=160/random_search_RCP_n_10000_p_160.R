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

source(here("lib", "simulation_routines.R"))
source(here("lib", "factorial_analysis_routines.R"))
source(here("lib", "gibbs_sampler_random_search.R"))


## Base directory to save MCMC results to
BASE_DIR = "simulation_factorial_analysis_n_10000_p_160_outputs"

## Read synthetic datasets for the 2 studies

results_n_10000_p_160 = readRDS(
  here("output",
       BASE_DIR,
       "results_n_10000_p_160.rds"))


X = results_n_10000_p_160.rds[['0.9']][['moderate']]$X.syn
y = results_n_10000_p_160.rds[['0.9']][['moderate']]$y.syn

## run the random search

v0 = rep(1, ncol(X))
chain <- suppressWarnings(
  Gibbs.random.search(v0, data.frame(X, y),
                      burnin=BURNIN,
                      iter=ITER,
                      penalty='BIC')
)


## check directory

cat(sprintf("[%s] Verifying existence of output directory...\n", Sys.time()))

outdir <- here("output", BASE_DIR)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat(sprintf("[%s] Output Directory Successfully Verified. \n", Sys.time()))



## save the chain
filename <- "gibbs_samples_n_10000_p_160.rds"
saveRDS(res.RS, file.path(outdir, filename))



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



## ad-hoc analysis
rel.freq = rel.freq.inference(chain)

MAP = find_MAP(chain)

BIC_minimizer = find_BIC_minimizer(data.frame(X,y), chain)













