#' Constructs a block diagonal correlation matrix used to generate
#' a synthetic multivariate normal design matrix later on.
#' For p predictors, it is assumed that X1 to X_{p/4} are ground
#' truth variables (involved in generating y) and X_{p/4 + 1 } to 
#' X_p are 'redundant' variables.
#' 
#' Out of the p/4 ground truths, the first half X_1 to X_{p/8} are correlated
#' (as in, mutual pairwise correlation)
#' and the second half X_{p/8 + 1} to X_{p/4} are not.
#'
#' @param p The number of predictors in total
#' @param rho Pairwise correlation coefficients for correlated ground truths
#'
#' @return A block diagonal matrix of size p with the above specifications.
construct_block_cov <- function(p,
                                rho) {
  

  p_corr <- 2
  p_uncorr_gt <- 2
  
  p_redundant <- p - p_corr - p_uncorr_gt
  
  # checks -- is p too small?
  if (p_corr + p_uncorr_gt > p) stop("p is too small to define ground truth blocks.")
  
  # correlated ground truths
  Sigma1 <- matrix(rho, nrow = p_corr, ncol = p_corr)
  diag(Sigma1) <- 1
  
  # generate the block diagonal matrix
  r <- as.matrix(Matrix::bdiag(Sigma1, diag(p_uncorr_gt), diag(p_redundant)))
  return(r)
}



# Construct a correlation matrix with GTs and correlated redundants
construct_corr_w_redundant <- function(p,
                                  gt_set     = 1:4,       # ground truth indices
                                  gt_corr_pairs = list(c(1,2)), # GTs to correlate
                                  rho_gt    = 0.9,        # corr inside GT pairs
                                  red_set   = 5:19,       # redundant indices
                                  gt_to_red = c(1,2),     # which GTs redundants tie to
                                  rho_gt_red= 0.9,        # corr GT↔redundants
                                  jitter_eps= 1e-8,
                                  nearest_PD=TRUE,
                                  variances=1) {
  stopifnot(p >= 1,
            all(gt_set >= 1 & gt_set <= p),
            all(red_set >= 1 & red_set <= p))
  
  Σ <- diag(p) * variances
  
  # 1) Correlate chosen GT pairs
  if (length(gt_corr_pairs)) {
    for (grp in gt_corr_pairs) {
      g <- intersect(gt_set, grp)
      if (length(g) >= 2) {
        for (i in 1:(length(g)-1)) for (j in (i+1):length(g)) {
          Σ[g[i], g[j]] <- Σ[g[j], g[i]] <- rho_gt
        }
      }
    }
  }
  
  # 2) Correlate redundants to specified GTs
  g2r <- intersect(gt_set, gt_to_red)
  rds <- intersect(red_set, 1:p)
  if (length(g2r) && length(rds)) {
    for (gi in g2r) for (rj in rds) {
      Σ[gi, rj] <- Σ[rj, gi] <- rho_gt_red
    }
  }
  
  # 3) note non zero determinant != positive definite (need all eigenvalues to be > 0)
  # if not PD, make an apprxoiamtion
  eig <- eigen(Σ, symmetric = TRUE, only.values = TRUE)$values
  cat("Smallest eigenvalue: ", min(eig), "\n")
  
  if (nearest_PD){
    if (min(eig) <= 1e-10) {
      if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
      if (!requireNamespace("sfsmisc", quietly = TRUE)) install.packages("sfsmisc")
      # nearest positive definite (alternative: Matrix::nearPD)
      Σ <- as.matrix(Matrix::nearPD(Σ, corr = TRUE)$mat)
      
      eig <- eigen(Σ, symmetric = TRUE, only.values = TRUE)$values
      cat("Smallest eigenvalue after trf: ", min(eig), "\n")
      
    }
  } else{
    # 3) Ensure symmetry + tiny jitter
    Σ <- (Σ + t(Σ)) / 2
    Σ <- Σ + diag(jitter_eps, p)
  }
  
  Σ
  
}


# How about we construct the covariance matrix
# by Σ = diag(s) * R * diag(s)
construct_corr_w_redundant_v2 <- function(
    p,
    gt_set        = 1:4,             # ground-truth indices
    gt_corr_pairs = list(c(1,2)),    # pairs/groups within gt_set to correlate
    rho_gt        = 0.9,             # correlation for GT↔GT pairs
    red_set       = 5:19,            # redundant indices
    gt_to_red     = c(1,2),          # which GTs redundants tie to
    rho_gt_red    = 0.9,             # correlation for GT↔redundant
    sds           = 1,               # scalar or length-p vector of SDs
    nearest_PD    = TRUE
) {
  stopifnot(p >= 1)
  # --- build correlation matrix R first ---
  R <- diag(p)
  # correlate specified GT pairs/groups
  if (length(gt_corr_pairs)) {
    for (grp in gt_corr_pairs) {
      g <- intersect(gt_set, grp)
      if (length(g) >= 2) {
        for (i in 1:(length(g)-1)) for (j in (i+1):length(g)) {
          R[g[i], g[j]] <- R[g[j], g[i]] <- rho_gt
        }
      }
    }
  }
  # correlate redundants to chosen GTs
  g2r <- intersect(gt_set, gt_to_red)
  rds <- intersect(red_set, 1:p)
  if (length(g2r) && length(rds)) {
    for (gi in g2r) for (rj in rds) {
      R[gi, rj] <- R[rj, gi] <- rho_gt_red
    }
  }
  
  # ensure positive definiteness (of R)
  eig_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (nearest_PD && eig_min <= 1e-10) {
    if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)  # stays a correlation matrix
  }
  
  # --- now scale to the desired variances without changing correlations ---
  if (length(sds) == 1L) sds <- rep(sds, p)
  D  <- diag(sds, p)
  Sigma <- D %*% R %*% D
  
  Sigma
}

# Build a PSD correlation (or covariance) matrix with
# strong GT↔redundant links AND induced within-redundant correlation.
construct_corr_w_redundant_v3 <- function(
    p,
    gt_set        = 1:4,                 # indices of ground truths
    gt_corr_pairs = list(c(1,2)),        # pairs inside gt_set to set at rho_gt
    rho_gt        = 0.9,                 # correlation for those pairs
    red_set       = 5:10,                # indices of redundants
    gt_to_red     = c(2,4),              # which GTs each redundant ties to
    rho_gt_red    = 0.9,                 # Corr(redundant, each chosen GT)
    eps_scale     = 1e-8,                # small slack for feasibility
    sds           = 1,                   # scalar or length-p vector; 1 => correlation matrix
    verbose       = TRUE
) {
  stopifnot(p >= 1, all(gt_set %in% seq_len(p)), all(red_set %in% seq_len(p)))
  gt_set <- sort(unique(gt_set))
  red_set <- sort(unique(red_set))
  other_set <- setdiff(seq_len(p), union(gt_set, red_set))
  
  g <- length(gt_set)
  m <- length(red_set)
  if (g == 0 || m == 0) stop("Need at least one GT and one redundant variable.")
  
  ## 1) GT–GT block G
  G <- diag(g)
  if (length(gt_corr_pairs)) {
    for (pair in gt_corr_pairs) {
      idx <- match(pair, gt_set)  # map global idx -> local GT idx
      idx <- idx[!is.na(idx)]
      if (length(idx) >= 2) {
        for (i in 1:(length(idx)-1)) for (j in (i+1):length(idx)) {
          G[idx[i], idx[j]] <- G[idx[j], idx[i]] <- rho_gt
        }
      }
    }
  }
  # quick PSD check / fix for G only (rarely needed if pairs are simple)
  evG <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  if (min(evG) <= 0) {
    if (verbose) message("Adjusting GT block to nearest PD.")
    G <- as.matrix(Matrix::nearPD(G, corr = TRUE)$mat)
  }
  
  ## 2) GT↔redundant links: rows of A are alpha_i (length g each)
  A <- matrix(0, nrow = m, ncol = g)
  tie_cols <- match(gt_to_red, gt_set)
  tie_cols <- tie_cols[!is.na(tie_cols)]
  if (length(tie_cols) == 0) stop("gt_to_red contains no indices inside gt_set.")
  A[, tie_cols] <- rho_gt_red
  
  ## 3) Enforce feasibility rowwise: s_i = alpha_i^T G^{-1} alpha_i ≤ 1 - eps
  Ginv <- solve(G)
  s_vec <- rowSums((A %*% Ginv) * A)  # diag(A G^{-1} A^T)
  over <- which(s_vec >= 1 - eps_scale)
  if (length(over)) {
    if (verbose) message("Scaling down ", length(over), " redundant link(s) to keep s_i <= 1.")
    for (i in over) {
      ci <- sqrt((1 - eps_scale) / s_vec[i])
      A[i, ] <- ci * A[i, ]
    }
    s_vec <- rowSums((A %*% Ginv) * A)  # recompute
  }
  
  ## 4) Redundant–redundant block: B = A G^{-1} A^T + diag(1 - s)
  B <- A %*% Ginv %*% t(A)
  diag(B) <- diag(B) + (1 - s_vec)  # now diag(B) = 1
  
  ## 5) Assemble full correlation matrix R
  R <- diag(p)
  # place G
  R[gt_set, gt_set] <- G
  # place A (GT↔red)
  R[red_set, gt_set] <- A
  R[gt_set, red_set] <- t(A)
  # place B
  R[red_set, red_set] <- B
  # everything else (others) stays independent
  
  # symmetry & final PSD check
  R <- (R + t(R)) / 2
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (verbose) message(sprintf("Smallest eigenvalue of correlation matrix: %.6f", min(ev)))
  
  ## 6) Return correlation or covariance (scale by sds)
  if (length(sds) == 1L) {
    if (sds == 1) return(R)  # correlation matrix
    D <- diag(sds, p)
    return(D %*% R %*% D)
  } else {
    stopifnot(length(sds) == p, all(sds > 0))
    D <- diag(sds, p)
    return(D %% R %*% D)
  }
}

# PSD correlation (or covariance) with GTs + correlated redundants
construct_corr_psd <- function(
    p,
    gt_set,                        # e.g. 1:4
    gt_corr_pairs = list(c(1,2)),  # pairs inside gt_set to set at rho_gt
    rho_gt = 0.9,
    red_set,                       # e.g. 5:15
    red_links,                     # list length = length(red_set);
    # each element: integer GT indices (global) or
    # a named/numeric vector c(gt_index = rho_gt_red)
    sds = 1,                       # 1 => correlation matrix; else covariance
    eps = 1e-8, verbose = TRUE
) {
  stopifnot(all(gt_set %in% 1:p), all(red_set %in% 1:p))
  gt_set <- sort(unique(gt_set)); red_set <- sort(unique(red_set))
  g <- length(gt_set); m <- length(red_set)
  if (length(red_links) != m) stop("red_links must match length(red_set).")
  
  # -- GT block G
  G <- diag(g)
  if (length(gt_corr_pairs)) {
    for (pair in gt_corr_pairs) {
      idx <- match(pair, gt_set); idx <- idx[!is.na(idx)]
      if (length(idx) >= 2)
        for (i in 1:(length(idx)-1)) for (j in (i+1):length(idx))
          G[idx[i], idx[j]] <- G[idx[j], idx[i]] <- rho_gt
    }
  }
  evG <- eigen(G, TRUE, TRUE)$values
  if (min(evG) <= 0) {
    if (verbose) message("Adjusting GT block to PD.")
    G <- as.matrix(Matrix::nearPD(G, corr = TRUE)$mat)
  }
  Ginv <- solve(G)
  
  # -- GT↔redundant matrix A (rows are alpha_i over GTs)
  A <- matrix(0, nrow = m, ncol = g)
  for (i in seq_len(m)) {
    li <- red_links[[i]]
    if (is.null(names(li))) {                 # only indices given
      cols <- match(li, gt_set)
      A[i, cols] <- 0.9                       # default strength
    } else {                                  # named/numeric vector: names are GT indices, values are rhos
      cols <- match(as.integer(names(li)), gt_set)
      A[i, cols] <- as.numeric(li)
    }
  }
  
  # -- Rowwise feasibility: alpha_i^T G^{-1} alpha_i <= 1 - eps
  s <- rowSums((A %*% Ginv) * A)
  over <- which(s >= 1 - eps)
  if (length(over)) {
    if (verbose) message("Scaling ", length(over), " redundant link(s) to keep feasibility.")
    for (i in over) A[i, ] <- A[i, ] * sqrt((1 - eps) / s[i])
    s <- rowSums((A %*% Ginv) * A)
  }
  
  # -- Redundant block B and full correlation
  B <- A %*% Ginv %*% t(A); diag(B) <- diag(B) + (1 - s)
  
  R <- diag(p)
  R[gt_set, gt_set]   <- G
  R[red_set, gt_set]  <- A
  R[gt_set, red_set]  <- t(A)
  R[red_set, red_set] <- B
  R <- (R + t(R)) / 2
  
  ev <- eigen(R, TRUE, TRUE)$values
  if (verbose) message(sprintf("min eigen(R) = %.6f", min(ev)))
  
  if (length(sds) == 1L && sds == 1) return(R)
  D <- diag(if (length(sds) == 1L) rep(sds, p) else sds, p)
  D %*% R %*% D
}




#' Simulate Binary Response and Assess Ground Truth Covariate Significance
#'
#' This function simulates a binary response variable `y` using a specified 
#' coefficient vector `beta` and a synthetic design matrix `X`. 
#' The first `num_ground_truths` covariates are assumed to be the true 
#' signal variables. A logistic model is fitted to these covariates and 
#' their significance is assessed.
#'
#' @param beta A numeric vector of coefficients used to generate the signal (true) covariates.
#' @param X A matrix or data frame of covariates (synthetic design matrix). ASSUMES columns 1 to `num_ground_truths` are ground truths.
#' @param num_ground_truths Number of ground truth covariates to use from `X`. Defaults to `length(beta)`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{X.syn}{The full synthetic design matrix (input).}
#'   \item{y.syn}{The simulated binary response vector.}
#'   \item{beta}{The coefficient vector used to generate `y`.}
#'   \item{mod}{The fitted logistic regression model (`glm` object).}
#'   \item{sig}{A data frame summarizing coefficient estimates, p-values, and significance of each ground truth variable.}
#' }
#'
#' @examples
#' r <- diag(4)
#' X <- MASS::mvrnorm(n = 2000, mu = rep(0, 4), Sigma = r, empirical = TRUE)
#' beta <- c(3, 2, -1.5, 2.5)
#' result <- check_significance(beta, X)
#'
#' @export
check_significance = function(beta,
                              X,
                              num_ground_truths=length(beta),
                              SEED=123,
                              plot=TRUE,
                              lower_prob=0.05,
                              upper_prob=0.95,
                              verbose=FALSE){
  
  ## 0. some helper functions
  
  signif_code <- function(p) {
    if (is.na(p)) return("NA")
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    if (p <= 0.1)   return(".")
    return(" ")
  }
  
  non_extreme <- function(probs, lower = lower_prob, upper = upper_prob) {
    mean(probs > lower & probs < upper)
  }
  
  set.seed(SEED)
  
  ## 1. generate responses based on the 
  ## given beta coefficients of GROUND truth variables
  ## by convention, denote X1 to X_{NUM_GROUND_TRUTHS} as the ground truths
  
  X.true = X[,1:num_ground_truths]
  
  #print(colnames(X.true))
  
  eta = X.true %*% beta
  p = 1 / (1 + exp(-eta))
  y = rbinom(nrow(X.true), size=1, prob=p)
  
  
  # see how significant each variable is by fitting a glm()
  if (verbose){
    cat("📍 Significance of Ground Truths Below: \n")
  }

  
  mod = glm(formula=y~.,
            family='binomial',
            data=data.frame(X.true, y))
  
  
  # summary of model is printed for your convenience
  if (verbose){
    print(summary(mod))
  }
  
  
  # track significance
  coef_summary <- summary(mod)$coefficients
  significance_table <- data.frame(
    Estimate = coef_summary[-1, "Estimate"],  # remove intercept
    P_Value = coef_summary[-1, "Pr(>|z|)"],
    Significant = coef_summary[-1, "Pr(>|z|)"]
  )
  
  ## 2. plot true probabilities
  probs_true <- plogis(eta)
  if(plot & verbose) {
    hist(probs_true, breaks=30, main="Histogram of True Probabilities", xlab="p_i = sigmoid(Xβ_true)")
    abline(v=c(0.1,0.9), col="red", lty=2)
  }
  
  ## 2.1 plot fitted proababilities
  betahat <- coef(mod)[-1]
  probs_fitted <- plogis(as.matrix(X.true) %*% betahat)
  if(plot & verbose) {
    hist(probs_fitted, breaks=30, main="Histogram of Fitted Probabilities", xlab="p_i = sigmoid(Xβ_hat)")
    abline(v=c(0.1,0.9), col="red", lty=2)
  }
  
  ## 2.2 Plot distribution of linear predictor
  if(plot & verbose) {
    hist(eta, breaks=30, main="Histogram of Linear Predictor", xlab="eta = Xβ_hat")
  }
  
  ## 2.2 Proportion of non-extreme probabilities
  prop_non_extreme_true <- non_extreme(probs_true)
  prop_non_extreme_fitted <- non_extreme(probs_fitted)
  
  if (verbose){
    cat(sprintf("📍 Fraction of true probabilities in (%.2f, %.2f): %.3f\n\n", 
                0.05, 0.95, prop_non_extreme_true))
    cat(sprintf("📍 Fraction of fitted probabilities in (%.2f, %.2f): %.3f\n\n", 
                0.05, 0.95, prop_non_extreme_fitted))    
  }
  

  
  
  ## 4. perform a LRT
  # Implementation: for each X1,..,X4, obtain its p-value 
  # from the LRT (not wald) test on H_0: beta_i = 0
  lrt_results <- sapply(1:num_ground_truths, function(j) {
    
    # obttain glm() formula string for full and reduced model
    # (1 variable less)
    vars <- paste0("X", 1:num_ground_truths)
    formula_full <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
    formula_reduced <- as.formula(paste("y ~", paste(vars[-j], collapse = " + ")))
    
    # fit both models
    mod_full <- glm(formula_full, family="binomial", data = data.frame(X.true, y))
    mod_reduced <- glm(formula_reduced, family="binomial", data = data.frame(X.true, y))
    
    # obtain p-value from LRT
    pval <- anova(mod_reduced, mod_full, test="LRT")[2, "Pr(>Chi)"]
    return(pval)
    
  })
  
  lrt_table <- data.frame(
    Variable = paste0("X", 1:num_ground_truths),
    LRT_P_Value = lrt_results
  )
  
  
  #cat("LRT P-values below:\n")
  #print(as.matrix(lrt_table))
  
  ## 4.1 Comparing LRT VS Wald P-values
  wald_pvals <- significance_table$P_Value
  wald_signif_codes <- sapply(wald_pvals, signif_code)
  lrt_pvals <- lrt_results
  lrt_signif_codes <- sapply(lrt_pvals, signif_code)
  
  lrt_vs_wald_table <- data.frame(
    Variable = paste0("X", 1:num_ground_truths),
    Wald_P = wald_pvals,
    wald_signif_codes = wald_signif_codes,
    LRT_P = lrt_pvals,
    lrt_signif_codes = lrt_signif_codes,
    LRT_minus_Wald = lrt_pvals - wald_pvals
  )
  
  lrt_vs_wald_table_2 <- data.frame(
    Variable = paste0("X", 1:num_ground_truths),
    LRT_minus_Wald = lrt_pvals - wald_pvals,
    LRT_div_Wald = lrt_pvals / wald_pvals,
    log10_LRT_div_Wald = log10(lrt_pvals / wald_pvals)
  )
  
  if (verbose){
    cat("📍 Difference b/w LRT VS Wald P-values below:\n\n")
    print(as.matrix(lrt_vs_wald_table))
    cat("\n")
    print(as.matrix(lrt_vs_wald_table_2))
    cat("\n")    
  }

  
  ## 5. final checks, return synthetic responses and their metadata
  
  # ensure distribution of y is not too skewed
  class_balance <- table(y) / length(y)
  
  if (verbose){
    cat("📍 Distribution of generated y: ", class_balance, " \n")
    
    if (any(class_balance < 0.2)) {
      warning("⚠️ Class imbalance detected: one class is <20% of total.")
    }
  }

  
  
  # output everything (including the input synthetic design matrix)
  out = list(
    X.syn=X,
    y.syn=y,
    beta=beta,
    mod=mod,
    sig=significance_table,
    lrt=lrt_table,
    lrt_vs_wald_table=lrt_vs_wald_table,
    probs_true=probs_true,
    probs_fitted=probs_fitted,
    class_balance=class_balance
  )
  
  return(out)
  
}


#' Simulates multivariate normal synthetic variables given a
#' covariance (or correlation) matrix r.
#' 
#' Main use case: generate candidate variables (ground truth + redundant)
#' for random search procedure. So dimension of r is the total number of 
#' candidate variables, NOT the number of ground truths.
#'
#'
#' @param r A pxp square, symmetric covariance matrix
#' @param SIZE The number of random vectors to be sampled
#'
#' @return Sampled random vectors stored in a SIZE x r matrix
generate_synthetic_vars = function(r, SIZE=2000, SEED=123){
  
  set.seed(SEED)
  # generates SIZE random vectors of dimension nrow(r) or ncol(r)
  syn.vars <- mvrnorm(n=SIZE,
                      mu=rep(0, nrow(r)),
                      Sigma=r,
                      empirical=TRUE )
  
  # include an intercept column? no need
  
  # fix column names
  varnames = 1:nrow(r)
  varnames = sapply(varnames, FUN=function(i){paste("X",i,sep='')})
  
  
  colnames(syn.vars) = varnames
  
  return(syn.vars)
  
}











