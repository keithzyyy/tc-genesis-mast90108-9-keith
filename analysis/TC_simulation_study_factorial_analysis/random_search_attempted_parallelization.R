# brainstorming as to how one can run all 15 simulation
# studies in parallel

# routine to run sim analysis. Requires cloning the following:
# "./output/simulation_factorial_analysis_outputs"
# source(here("lib", "simulation_routines.R"))
# source(here("lib", "factorial_analysis_routines.R"))
# source(here("lib", "gibbs_sampler_random_search.R"))


run_simulation_analysis <- function(i) {
  library(here)
  
  # Read in the result
  result <- readRDS(
    here("output", "simulation_factorial_analysis_outputs", paste0("result.c", i, ".rds"))
  )
  
  X <- result$X.syn
  y <- result$y.syn
  
  suppressWarnings({
    gibbs.RS <- Gibbs.random.search.basin(
      rep(1, ncol(X)),
      data.frame(cbind(X, y)),
      burnin = 200,
      iter = 1000,
      penalty = 'BIC'
    )
  })
  
  rel.freqs <- rel.freq.inference(gibbs.RS, CUTOFF = 0.5)
  best.Xs <- colnames(X)[as.logical(find_MAP(gibbs.RS))]
  
  # Save outputs
  saveRDS(gibbs.RS, here("output", "simulation_factorial_analysis_outputs", paste0("gibbs.RS.c", i, ".rds")))
  saveRDS(rel.freqs, here("output", "simulation_factorial_analysis_outputs", paste0("gibbs.RS.c", i, ".rel.freqs.rds")))
  saveRDS(best.Xs, here("output", "simulation_factorial_analysis_outputs", paste0("gibbs.RS.c", i, ".best.Xs.rds")))
  
  return(paste("Completed i =", i))
}


library(future.apply)
plan(multisession, workers = 4)  # or multicore if on Unix

results <- future_lapply(1:15, run_simulation_analysis)
