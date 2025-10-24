# Load all local and ext libraries ----
if (!requireNamespace("here", quietly=TRUE)) install.packages("here")
library(here) # to override here(), simply go to File -> New Project
here()
library(parallel)

# one hot encoding with intercept
source(here("lib", "preprocessing_v3.R"))
# full rank encoding with NO intercept
source(here("lib", "preprocessing_v2.R"))
source(here("lib", "simulation_routines.R"))
source(here("lib", "factorial_analysis_routines.R"))
source(here("lib", "gibbs_sampler_random_search.R"))

# Essential params ----
SEED=123
ITER   <- 30200
BURNIN <- 200
GLOBAL_START_TIME <- Sys.time()
set.seed(SEED) # NOTE: NO further set.seed() needs to be invoked.
NUM_CHAINS <- 4


# Prep the TC genesis dataset ----
X_TC_FP = here("data", "curated","X_TC_genesis.rds")
Y_TC_FP = here("data", "curated","y_TC_genesis.rds")
X_TC_EXIST = file.exists(X_TC_FP)
Y_TC_EXIST = file.exists(Y_TC_FP)

if (X_TC_EXIST & Y_TC_EXIST){
  cat("Preprocessed (X,y) already exist. Importing now.\n")
  X = readRDS(X_TC_FP)
  y = readRDS(Y_TC_FP)
} else {
  cat("Preprocessing the TC genesis dataset....\n")
    data_path = here("data", "envDataset_12h_10x10_with_mask.csv")
    # import the data
    data = read.csv(data_path)
    # fix basin NA category first (NA missing to "NA")
    data.NA.fixed = fix_NA_value_and_refactor_basin(data)
    X = data.NA.fixed[, ! colnames(data.NA.fixed) %in% ( c("basin", "TC_genesis") )]
    y = data.NA.fixed$TC_genesis
    
    saveRDS(X, here("data", "curated","X_TC_genesis.rds"))
    saveRDS(y, here("data", "curated","y_TC_genesis.rds"))
}

# Functions: tuning lambda in MH----

tune_lambda <- function(lambdas, v0, data, steps = 500, burnin = 100, seed = 123) {
  set.seed(seed)
  stopifnot(length(lambdas) > 0, steps > burnin)
  
  start_time = Sys.time()
  
  
  out <- lapply(lambdas, function(lam) {
    cat(sprintf("λ = %s ...\n", lam))
    rs <- MH.random.search(v0, data, iter = steps, burnin = burnin,
                           lambda = lam, debug = FALSE)
    samp <- rs$samples
    # distinct-model proportion
    keys <- apply(samp, 1, paste0, collapse = "")
    distinct_prop <- length(unique(keys)) / nrow(samp)
    
    cat("Acceptance rate ", rs$acc.rate, "\n")
    cat("Distinct models ", distinct_prop, "\n")
    
    data.frame(lambda = lam, acc = rs$acc.rate, distinct = distinct_prop)
  })
  
  cat("Duration (min):", round(as.numeric(difftime(Sys.time(), start_time, units="mins")), 2), "\n")
  do.call(rbind, out)
}


# Parallel λ–tuning with MH.random.search + progress bar

tune_lambda_parallel <- function(lambdas, v0, data,
                                 steps = 500, burnin = 100,
                                 seed = 123, workers = NULL) {
  # Basic sanity checks
  stopifnot(length(lambdas) > 0, steps > burnin)
  
  cat(
    sprintf("Tuning temp for RS MH, each with %d models excluding %d burnin\n",
            steps, burnin)
    )
  
  # Ensure needed packages are installed
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
  if (!requireNamespace("progressr", quietly = TRUE))    install.packages("progressr")
  
  # Generate a reproducible seed for each λ
  lam_seeds <- seed + seq_along(lambdas) * 1000L
  
  # Save the current parallel plan, so we can restore it at the end
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Default number of workers = min(#λ values, available CPU cores)
  if (is.null(workers)) {
    workers <- min(length(lambdas), parallel::detectCores(logical = TRUE))
    cat(sprintf("%d workers used for %d lambda values\n", workers, length(lambdas)))
  }
  # Launch multisession workers (works on Windows, macOS, Linux)
  future::plan(future::multisession, workers = workers)
  
  # Set up progress bar
  progressr::handlers(global = TRUE)   # choose your favorite handler style globally
  
  res <- NULL
  progressr::with_progress({
    # Initialize progress bar with #λ steps
    p <- progressr::progressor(steps = length(lambdas))
    
    t0 <- Sys.time()
    # Loop in parallel over λ values
    res <- future.apply::future_mapply(
      FUN = function(lam, sd) {
        # Reproducibility per λ
        set.seed(sd)
        
        # Run MH random search for this λ
        rs <- MH.random.search(v0, data, iter = steps, burnin = burnin,
                               lambda = lam, debug = FALSE)
        
        # Proportion of distinct models visited
        keys <- apply(rs$samples, 1, paste0, collapse = "")
        distinct_prop <- length(unique(keys)) / nrow(rs$samples)
        
        # Update progress bar with inline message
        p(sprintf("λ = %.3g | acc = %.3f | distinct = %.3f",
                  lam, rs$acc.rate, distinct_prop))
        
        # Return tidy row for this λ
        data.frame(lambda   = lam,
                   acc      = rs$acc.rate,
                   distinct = distinct_prop)
      },
      lam = lambdas, sd = lam_seeds,
      SIMPLIFY = FALSE
    )
    # Print runtime in minutes
    cat("Duration (min):",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
  })
  
  # Bind all rows into a single data.frame
  do.call(rbind, res)
}


# --- inputs --------------------------------------------------------------
v0 <- rep(1, ncol(X))

# numeric vector of λ to try
lambdas <- seq(0.5, 11, by = 0.25)

# --- run --------------------------------------------------------------
MH_tune_lambda <- tune_lambda_parallel(
  lambdas, v0, data.frame(X, y),
  steps = 2300, burnin = 300,
  seed = SEED
)

# --- save --------------------------------------------------------------
outdir <- here("output", "random_search_TC_dataset_tuning_outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
saveRDS(MH_tune_lambda, file = file.path(outdir, "mh_lambda_tuning_results.rds"))


#head(MH_tune_lambda[order(MH_tune_lambda$distinct, decreasing = TRUE),])

MIN_N_DISTINCT = 0.05

MH_tune_lambda_best_lam = MH_tune_lambda %>% filter(distinct >= MIN_N_DISTINCT)

# plot tuning results ----
op <- par(mfrow = c(1,2))

plot(MH_tune_lambda$lambda, MH_tune_lambda$acc, type="b",
     xlab=expression(lambda), ylab="Acceptance rate",
     xlim=c(0,11))

plot(MH_tune_lambda$lambda, MH_tune_lambda$distinct, type="b",
     xlab=expression(lambda), ylab="Distinct models (proportion)",
     xlim=c(0,11))

lines(MH_tune_lambda_best_lam$lambda, MH_tune_lambda_best_lam$distinct, type="b",
     xlab=expression(lambda), ylab="Distinct models (proportion)",
     xlim=c(0,11),
     col='red')

par(op)


# if we wanna choose lambda s.t. number of distinct models >= 0.05*K,
# then lambda = either 0.5 or 0.75.
# I'd say pretty reasonable, since we wanna explore. 


