# Load all local and ext libraries ----
if (!requireNamespace("here", quietly=TRUE)) install.packages("here")
library(here) # to override here(), simply go to File -> New Project
here()
library(parallel)
library(future.apply)
library(progressr)

source(here("lib", "gibbs_sampler_random_search.R"))
source(here("lib", "load_and_preprocess_TC.R"))

# Essential params ----
SEED=123
ITER   <- 41000
BURNIN <- 1000
GLOBAL_START_TIME <- Sys.time()
set.seed(SEED) # NOTE: NO further set.seed() needs to be invoked.



# 1. build diverse initial states (same across λ) ----
# make init states as overdispersed as possible 
# e.g. vec of all 1s, all 0s, sparse, random
make_inits <- function(p) {
  list(
    full      = rep(1L, p),
    empty     = rep(0L, p),
    sparse  = rbinom(p, 1, 0.05),
    random  = rbinom(p, 1, 0.5)
  )
}

# 2. grid: reuse the same 4 seeds for all λ ----
# 1. specify the 4 initial values
lambdas <- c(1.5, 1.75, 2.5, 3)
chains_per_lambda <- 4

# e.g. if SEED=123 then 124, 125, 126, 127
base_chain_seeds <- SEED + seq_len(chains_per_lambda)      # same across λ
p <- ncol(X)

init_list <- make_inits(p)                                      # same across λ
stopifnot(length(init_list) == chains_per_lambda)

cat("Running multiple RS with lambdas=",lambdas,"\n")

grid <- do.call(rbind, lapply(lambdas, function(lam) {
  data.frame(lambda = lam,
             chain  = seq_len(chains_per_lambda),
             seed   = base_chain_seeds,
             init_id= names(init_list),
             stringsAsFactors = FALSE)
}))

# for example
#> grid
#lambda chain seed init_id
#1    0.50     1  124    full
#2    0.50     2  125   empty
#3    0.50     3  126  sparse
#4    0.50     4  127  random
#5    0.75     1  124    full
#6    0.75     2  125   empty
#7    0.75     3  126  sparse
#8    0.75     4  127  random
#9    1.00     1  124    full
#10   1.00     2  125   empty
#11   1.00     3  126  sparse
#12   1.00     4  127  random
#13   1.25     1  124    full
#14   1.25     2  125   empty
#15   1.25     3  126  sparse
#16   1.25     4  127  random



# 3. Function to run ONE MH ----
run_MH_chain <- function(v0, data, iter, burnin, lambda, seed) {
  set.seed(seed)
  MH.random.search(v0, data, iter = iter, burnin = burnin,
                   lambda = lambda, debug = FALSE)
}



# 4. parallel execution with progress --------------------------------------

# build the parallel setup 
old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
workers <- min(nrow(grid), parallel::detectCores(TRUE))
future::plan(future::multisession, workers = workers)

progressr::handlers(global = TRUE)

outdir <- here::here("output","mh_runs_4_lambda_run_2")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

progressr::with_progress({
  pbar <- progressr::progressor(steps = nrow(grid))
  
  results <- future.apply::future_lapply(seq_len(nrow(grid)), function(i) {
    g  <- grid[i, ]
    v0 <- init_list[[ g$init_id ]]
    
    rs <- run_MH_chain(v0, data.frame(X, y),
                       iter = ITER, burnin = BURNIN,
                       lambda = g$lambda, seed = g$seed)
    
    pbar(sprintf("λ=%.2f | chain=%d (%s) | acc=%.3f",
                 g$lambda, g$chain, g$init_id, rs$acc.rate))
    
    # ---- summary for quick scanning
    distinct_prop <- {
      k <- apply(rs$samples, 1, paste0, collapse = "")
      length(unique(k)) / nrow(rs$samples)
    }
    mean_active <- mean(rowSums(rs$samples))
    summ <- list(
      lambda = g$lambda,
      chain  = g$chain,
      init   = g$init_id,
      acc    = rs$acc.rate,
      distinct_prop = distinct_prop,
      mean_active   = mean_active
    )
    
    # ---- save once, unified structure
    fn <- sprintf("mh_lambda_%0.2f_chain_%d_%s.rds",
                  g$lambda, g$chain, g$init_id)
    
    # save the entire rs output since it also contains the computed BICs
    saveRDS(list(meta = g, summary = summ, rs=rs),
            file.path(outdir, fn))
    
    TRUE
  }, future.seed = TRUE)
})



## Save metadata for reproducibility ----
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



