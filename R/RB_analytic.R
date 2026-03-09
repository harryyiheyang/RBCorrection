# ============================================================
# Analytic Rao-Blackwellization for onlyexposure = TRUE (p <= ~10)
# Uses df-shift derivative identities:
#   P'(kappa)  = 0.5 * (P_{df+2}(kappa) - P_{df}(kappa))
#   P''(kappa) = 0.25 * (P_{df+4} - 2 P_{df+2} + P_df)
# so g1 = d/dkappa log P, g2 = d^2/dkappa^2 log P
# ============================================================

matrixsqrt=function(A){
  fit=matrixEigen(A)
  d=c(fit$value)
  d1=d*0
  d1[d>0]=1/d[d>0]
  d=sqrt(d)
  d1=sqrt(d1)
  A=matrixMultiply(fit$vector,t(fit$vector)*d)
  B=matrixMultiply(fit$vector,t(fit$vector)*d1)
  C=list(w=A,wi=B,eigenfit=fit)
  return(C)
}
# ---- compute g1, g2 from noncentral chi-square tails via pchisq ----
g1g2_from_pchisq <- function(df, kappa, t) {
  stopifnot(df >= 1, kappa >= 0, t >= 0)

  # log upper-tail probabilities for stability
  logP0 <- pchisq(t, df = df,     ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP2 <- pchisq(t, df = df + 2, ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP4 <- pchisq(t, df = df + 4, ncp = kappa, lower.tail = FALSE, log.p = TRUE)

  if (!is.finite(logP0)) {
    return(list(P = 0, logP = logP0, g1 = NA_real_, g2 = NA_real_,
                r2 = NA_real_, r4 = NA_real_))
  }

  # ratios P_{df+2}/P_df and P_{df+4}/P_df computed in log-scale
  r2 <- exp(logP2 - logP0)
  r4 <- exp(logP4 - logP0)

  # g1 = (P'/P) = 0.5*(P2/P0 - 1)
  g1 <- 0.5 * (r2 - 1)

  # P''/P = 0.25*(P4/P0 - 2*P2/P0 + 1)
  P2_over_P <- 0.25 * (r4 - 2 * r2 + 1)

  # g2 = (P''/P) - (P'/P)^2
  g2 <- P2_over_P - g1^2

  list(P = exp(logP0), logP = logP0, g1 = g1, g2 = g2, r2 = r2, r4 = r4)
}

# ---- one-row analytic RB (matches your C++ mapping) ----
rb_onlyexposure_analytic_one <- function(beta_hat, Sigma_full, eta, cutoff, p = length(beta_hat) - 1) {
  beta_hat <- as.numeric(beta_hat)
  p1 <- length(beta_hat)
  stopifnot(is.matrix(Sigma_full), nrow(Sigma_full) == p1, ncol(Sigma_full) == p1)
  stopifnot(p >= 1, p < p1, eta > 0, cutoff >= 0)

  bx  <- beta_hat[1:p]
  Sxx <- Sigma_full[1:p, 1:p, drop = FALSE]

  # stable inverse of Sxx
  cholS <- chol(Sxx)
  Sxx_inv <- chol2inv(cholS)

  # kappa and t for noncentral chi-square tail
  q     <- as.numeric(t(bx) %*% Sxx_inv %*% bx)
  kappa <- q / (eta^2)
  tval  <- cutoff / (eta^2)

  # g1, g2
  gg <- g1g2_from_pchisq(df = p, kappa = kappa, t = tval)
  Pacc <- gg$P
  if (!is.finite(Pacc) || Pacc <= 0) {
    stop("Selection probability under the model is 0 or non-finite. Check cutoff/eta/Sigma/beta.")
  }
  g1 <- gg$g1
  g2 <- gg$g2

  # Conditional moments for exposure-part noise e_x | A, beta:
  # E(e_x|A)   = 2 g1 * beta_x
  # Cov(e_x|A) = eta^2(1+2g1) Sxx + 4 g2 beta_x beta_x'
  Ex   <- 2 * g1 * bx
  Covx <- eta^2 * (1 + 2 * g1) * Sxx + 4 * g2 * tcrossprod(bx)

  # Lift to full vector e = (e_x, e_r) using Gaussian conditioning
  Sxy <- Sigma_full[1:p, (p+1):p1, drop = FALSE]
  Syx <- Sigma_full[(p+1):p1, 1:p, drop = FALSE]
  Syy <- Sigma_full[(p+1):p1, (p+1):p1, drop = FALSE]

  Amap <- Syx %*% Sxx_inv
  Er   <- as.numeric(Amap %*% Ex)

  Schur <- Syy - Amap %*% Sxy
  Covr  <- Amap %*% Covx %*% t(Amap) + eta^2 * Schur
  Covxr <- Covx %*% Sxx_inv %*% Sxy

  # Assemble mean/cov of e | A
  E_e <- c(Ex, Er)

  Cov_e <- matrix(0, p1, p1)
  Cov_e[1:p, 1:p] <- Covx
  Cov_e[(p+1):p1, (p+1):p1] <- Covr
  Cov_e[1:p, (p+1):p1] <- Covxr
  Cov_e[(p+1):p1, 1:p] <- t(Covxr)

  # Map to init estimator: beta_init = beta_hat - e/eta^2
  beta_RB <- beta_hat - E_e / (eta^2)

  # Conditional covariance of beta_init given selection (what your C++ cov(acc_mat) estimates)
  cov_init_cond <- Cov_e / (eta^4)

  # Final RB covariance you compute in R:
  # Sigma_RB = (1 + 1/eta^2) Sigma_full - cov_init_cond
  cov_RB <- (1 + 1/eta^2) * Sigma_full - cov_init_cond

  # Symmetrize for numerical stability
  cov_init_cond <- 0.5 * (cov_init_cond + t(cov_init_cond))
  cov_RB        <- 0.5 * (cov_RB + t(cov_RB))

  list(
    beta_RB       = beta_RB,
    cov_init_cond = cov_init_cond,
    cov_RB        = cov_RB,
    se_RB         = sqrt(pmax(diag(cov_RB), 0)),
    P_accept      = Pacc,
    kappa         = kappa,
    t             = tval,
    g1            = g1,
    g2            = g2
  )
}

# ---- batch wrapper (m rows) with optional fallback logic like your min_accept ----
rb_onlyexposure_analytic_batch <- function(beta_select, SE_Select, Rxy, ldsc = rep(1,nrow(beta_select)), gcov = matrix(0,nrow=ncol(beta_select),ncol=ncol(beta_select)),
                                           eta, pv.threshold,kappa_thres = 3) {
  beta_select <- as.matrix(beta_select)
  m  <- nrow(beta_select)
  p1 <- ncol(beta_select)
  p = ncol(beta_select) - 1
  cutoff = qchisq(pv.threshold, p, lower.tail = FALSE)

  SE_Select <- as.matrix(SE_Select)
  Rxy <- as.matrix(Rxy)
  RxyList=RxysqrtList=array(0,c(p1,p1,m))
  for(i in 1:m){
    sei=SE_Select[i,]
    Rxyi=t(t(Rxy)*sei)*sei+ldsc[i]*gcov
    RxyList[,,i]=Rxyi
    RxysqrtList[,,i]=matrixsqrt(Rxyi)$w
  }

  stopifnot(is.array(RxyList), length(dim(RxyList)) == 3,
            dim(RxyList)[1] == p1, dim(RxyList)[2] == p1, dim(RxyList)[3] == m)

  BETA_RB <- beta_select
  SE_RB   <- matrix(NA_real_, m, p1)
  COV_RB  <- vector("list", m)
  corrected <- integer(0)

  for (i in seq_len(m)) {
    Sigma_full <- RxyList[, , i]

    out <- rb_onlyexposure_analytic_one(beta_hat = beta_select[i, ],
                                        Sigma_full = Sigma_full,
                                        eta = eta, cutoff = cutoff, p = p)

    # Mimic C++: max_draws = 5*B proposals, require >= min_accept accepted
    # Expected accepted ≈ 5*B*P_accept
    if (is.finite(out$P_accept) && (out$P_accept >= 0.1) ) {
      BETA_RB[i, ] <- out$beta_RB
      SE_RB[i, ]   <- out$se_RB
      COV_RB[[i]]  <- out$cov_RB
      corrected <- c(corrected, i)
    } else {
      # fallback (same spirit as your code)
      BETA_RB[i, ] <- beta_select[i, ]
      SE_RB[i, ]   <- NA_real_
      M <- matrix(0, p1, p1)
      M[p1, p1] <- -1
      COV_RB[[i]] <- M
    }
  }


  list(
    BETA_RB = BETA_RB,
    bX_RB = BETA_RB[corrected,1:p],
    by_RB = BETA_RB[corrected,p+1],
    bXse_RB = SE_RB[corrected,1:p],
    byse_RB = SE_RB[corrected,p+1],
    COV_RB = COV_RB[corrected],
    CORRECTED_INDICES = corrected
  )
}

# ============================================================
# Example usage:
# res <- rb_onlyexposure_analytic_batch(beta_select, RxyList, eta, cutoff, B=1000, min_accept=100)
# ============================================================
