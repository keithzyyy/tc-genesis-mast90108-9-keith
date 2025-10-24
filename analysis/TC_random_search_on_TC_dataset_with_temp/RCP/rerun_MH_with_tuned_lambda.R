# Load all local and ext libraries ----
if (!requireNamespace("here", quietly=TRUE)) install.packages("here")
library(here) # to override here(), simply go to File -> New Project
here()
library(parallel)


source(here("lib", "gibbs_sampler_random_search.R"))
source(here("lib", "load_and_preprocess_TC"))

# Essential params ----
SEED=123
ITER   <- 30200
BURNIN <- 200
GLOBAL_START_TIME <- Sys.time()
set.seed(SEED) # NOTE: NO further set.seed() needs to be invoked.
NUM_CHAINS <- 4