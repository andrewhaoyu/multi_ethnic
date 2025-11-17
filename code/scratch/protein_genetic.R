# Minimal, runnable demo: joint G + X with partial observability (two-stage only)
# ------------------------------------------------------------
# What this script does (simple and stable):
# 1) Simulates:
#    - Genotypes G (light LD via AR(1) correlation and Gaussian copula → {0,1,2}).
#    - Proteins X | G = G %*% Gamma + MVN(0, Sigma_X).
#    - Outcome Y = G %*% beta_G + X %*% beta_X + e  (Gaussian) OR logistic for binary.
#    - Only n2 subjects have X observed; others are missing.
# 2) Two-stage estimator:
#    - Stage 1: Fit X ~ G with multi-response lasso (glmnet mgaussian) on observed-X subset; impute Xhat for all.
#    - Stage 2: Fit Y ~ [G, X_mixed] where X_mixed uses observed X when available else Xhat.
# 3) Prints test MSE (Gaussian) or AUC (Binary) and prevalence.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(MASS)      # mvrnorm
  library(glmnet)    # lasso / elastic net
  library(Matrix)
  library(pROC)      # AUC
})

# -----------------------------
# Helpers: LD correlation and genotype simulator
# -----------------------------
make_block_ar1_corr <- function(p, n_blocks = 5, rho = 0.5) {
  stopifnot(p >= n_blocks)
  base <- floor(p / n_blocks)
  sizes <- rep(base, n_blocks)
  rmd <- p - sum(sizes)
  if (rmd > 0) sizes[seq_len(rmd)] <- sizes[seq_len(rmd)] + 1
  blocks <- lapply(sizes, function(m){
    i <- seq_len(m); j <- seq_len(m)
    outer(i, j, function(a,b) rho^abs(a-b))
  })
  R <- as.matrix(Matrix::bdiag(blocks))
  diag(R) <- 1
  R
}

simulate_genotypes <- function(n, p, R, maf_low = 0.05, maf_high = 0.45) {
  stopifnot(is.matrix(R) && nrow(R) == p && ncol(R) == p)
  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = R)
  mafs <- runif(p, maf_low, maf_high)
  cut0 <- qnorm((1 - mafs)^2)
  cut1 <- qnorm(1 - mafs^2)
  G <- matrix(0, n, p)
  for (k in seq_len(p)) G[,k] <- (Z[,k] > cut0[k]) + (Z[,k] > cut1[k])
  storage.mode(G) <- "double"
  colnames(G) <- paste0("SNP", seq_len(p))
  G
}

# -----------------------------
# Data generation (small, simple)
# -----------------------------
simulate_small <- function(n1 = 500, n2 = 150, p = 30, q = 6,
                           outcome = c("gaussian","binomial"),
                           beta0 = -1.8, ld_blocks = 5, ld_rho = 0.5,
                           sigma2_X = 1, rho_X = 0.3, sigma2_e = 1, seed = 1) {
  outcome <- match.arg(outcome)
  set.seed(seed)
  
  # G
  R <- make_block_ar1_corr(p, n_blocks = ld_blocks, rho = ld_rho)
  G <- simulate_genotypes(n1, p, R)
  
  # Parameters (sparse but simple)
  set.seed(seed + 1)
  Gamma <- matrix(0, p, q)
  for (j in seq_len(q)) Gamma[sample.int(p, 6), j] <- rnorm(6, 0, 0.10)
  beta_G <- rep(0, p); beta_G[sample.int(p, 8)] <- rnorm(8, 0, 0.06)
  beta_X <- rep(0, q); beta_X[sample.int(q, 3)] <- rnorm(3, 0, 0.15)
  
  # Sigma_X and X
  i <- seq_len(q); j <- seq_len(q)
  Sigma_X <- sigma2_X * rho_X^abs(outer(i,j,"-")); diag(Sigma_X) <- sigma2_X
  X_full <- G %*% Gamma + MASS::mvrnorm(n1, rep(0, q), Sigma_X)
  colnames(X_full) <- paste0("Prot", seq_len(q))
  
  # Outcome
  lin <- as.vector(G %*% beta_G + X_full %*% beta_X)
  if (outcome == "gaussian") {
    y <- lin + rnorm(n1, 0, sqrt(sigma2_e))
  } else {
    prob <- plogis(beta0 + lin)
    y <- rbinom(n1, 1, prob)
  }
  
  # Partial X
  idx_obs <- sort(sample.int(n1, n2, replace = FALSE))
  X_obs <- matrix(NA_real_, n1, q); X_obs[idx_obs, ] <- X_full[idx_obs, ]
  
  list(G = G, X_obs = X_obs, X_full = X_full, y = y, idx_obs = idx_obs,
       outcome = outcome)
}

# -----------------------------
# Two-stage estimator
# -----------------------------
fit_stage1_X_given_G <- function(G, X_obs) {
  obs <- which(rowSums(is.finite(X_obs)) == ncol(X_obs))
  stopifnot(length(obs) >= 10)
  cvfit <- cv.glmnet(G[obs, , drop=FALSE], X_obs[obs, , drop=FALSE],
                     family = "mgaussian", alpha = 1, nfolds = 5)
  predict_X <- function(newG) {
    out <- drop(predict(cvfit, newx = newG, s = "lambda.min"))
    if (is.null(dim(out))) out <- matrix(out, nrow(newG))
    out
  }
  list(cvfit = cvfit, predict_X = predict_X, obs_index = obs)
}

make_X_mixed <- function(X_obs, Xhat) {
  Xmix <- Xhat
  obs <- which(rowSums(is.finite(X_obs)) == ncol(X_obs))
  if (length(obs) > 0) Xmix[obs, ] <- X_obs[obs, ]
  Xmix
}

fit_stage2_Y_on_G_X <- function(G, Xmix, y, family) {
  fam <- if (family == "gaussian") "gaussian" else "binomial"
  cvfit <- cv.glmnet(cbind(G, Xmix), y, family = fam, alpha = 1, nfolds = 5)
  list(cvfit = cvfit)
}

predict_stage2 <- function(fit, G, Xmix, family) {
  type <- if (family == "gaussian") "response" else "response"
  drop(predict(fit$cvfit, newx = cbind(G, Xmix), s = "lambda.min", type = type))
}

# -----------------------------
# One-click demo
# -----------------------------
run_simple_demo <- function(outcome = c("gaussian","binomial"), seed = 42) {
  outcome <- match.arg(outcome)
  cat("\n=== Running simple two-stage demo (", outcome, ") ===\n", sep = "")
  dat <- simulate_small(outcome = outcome, seed = seed)
  
  # Train/test split
  n1 <- nrow(dat$G)
  set.seed(seed + 100)
  tr <- sort(sample.int(n1, floor(0.7 * n1)))
  te <- setdiff(seq_len(n1), tr)
  
  Gtr <- dat$G[tr, , drop=FALSE]; Gte <- dat$G[te, , drop=FALSE]
  Xobstr <- dat$X_obs[tr, , drop=FALSE]; Xobste <- dat$X_obs[te, , drop=FALSE]
  ytr <- dat$y[tr]; yte <- dat$y[te]
  
  # Stage 1
  s1 <- fit_stage1_X_given_G(Gtr, Xobstr)
  Xhat_tr <- s1$predict_X(Gtr)
  Xhat_te <- s1$predict_X(Gte)
  
  # Build mixed X
  Xmix_tr <- make_X_mixed(Xobstr, Xhat_tr)
  Xmix_te <- make_X_mixed(Xobste, Xhat_te)
  
  # Stage 2
  s2 <- fit_stage2_Y_on_G_X(Gtr, Xmix_tr, ytr, family = outcome)
  pred <- predict_stage2(s2, Gte, Xmix_te, family = outcome)
  
  if (outcome == "gaussian") {
    mse <- mean((yte - pred)^2)
    cat("Test MSE:", sprintf("%.4f", mse), "\n")
  } else {
    auc <- as.numeric(pROC::auc(yte, pred))
    cat("Test AUC:", sprintf("%.4f", auc), "\n")
    cat("Prevalence (test):", sprintf("%.3f", mean(yte)), "\n")
  }
  
  invisible(list(stage1 = s1, stage2 = s2))
}

# -----------------------------
# Example runs
# -----------------------------
if (sys.nframe() == 0) {
  set.seed(123)
  run_simple_demo("gaussian", seed = 11)
  run_simple_demo("binomial", seed = 21)
}
