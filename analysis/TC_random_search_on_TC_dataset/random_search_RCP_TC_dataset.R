#!/usr/bin/env Rscript

## set metadata 
SEED=123
ITER   <- 2300
BURNIN <- 300
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


## Load all local libraries

# one hot encoding with intercept
source(here("lib", "preprocessing_v3.R"))
# full rank encoding with NO intercept
source(here("lib", "preprocessing_v2.R"))
source(here("lib", "simulation_routines.R"))
source(here("lib", "factorial_analysis_routines.R"))
source(here("lib", "gibbs_sampler_random_search.R"))

## Prep the TC genesis dataset

data_path = here("data", "envDataset_12h_10x10_with_mask.csv")
# import the data
data = read.csv(data_path)
# fix basin NA category first (NA missing to "NA")
data.NA.fixed = fix_NA_value_and_refactor_basin(data)


# note the only categorical variable is `basin` and we wanna include
# all basin levels in the model
X = data.NA.fixed[, ! colnames(data.NA.fixed) %in% ( c("basin", "TC_genesis") )]
y = data.NA.fixed$TC_genesis


## Run the random search
v0 = rep(1, ncol(X))
cat(sprintf("Running on %d rows and %d predictors.\n", nrow(X), ncol(X)))
chain = Gibbs.random.search(v0,
                              data.frame(X, y),
                              burnin=BURNIN,
                              iter=ITER,
                              penalty='BIC')

## Verify existence of output directory
cat(sprintf("[%s] Verifying existence of output directory...\n", Sys.time()))
BASE_DIR = "random_search_TC_dataset_outputs"
outdir <- here("output", BASE_DIR)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
cat(sprintf("[%s] Output Directory Successfully Verified. \n", Sys.time()))


## Save the mcmc samples
filename <- "gibbs_TC_dataset_results.rds"
saveRDS(chain, file.path(outdir, filename))
cat(sprintf("[%s] Saved %s\n", Sys.time(), filename))


## Save metadata for reproducibility
GLOBAL_END_TIME = Sys.time()
GLOBAL_TIME_TAKEN = as.numeric(difftime(GLOBAL_END_TIME,
                                        GLOBAL_START_TIME, units="mins"))

metadata <- list(
  SEED             = SEED,
  ITER             = ITER,
  BURNIN           = BURNIN,
  n_rows = nrow(X),
  n_cols = ncol(X),
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


