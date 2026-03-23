#' Rao-Blackwell Correction for Selected IVs in MVMR
#'
#' Performs Rao-Blackwell correction for instrumental variables using an analytic
#' closed-form solution based on noncentral chi-square distribution theory.
#'
#' @param BETA_Select A numeric matrix of dimension m x p+1: effect estimates.
#' @param SE_Select A numeric matrix of same dimension: standard errors.
#' @param Rxy A (p+1) x (p+1) correlation matrix of Z-statistics.
#' @param gcov A (p+1) x (p+1) matrix of per-SNP genetic covariances.
#' @param ldsc A numeric vector of length m: LD scores for each SNP.
#' @param eta Randomization parameter controlling noise level (default = 0.5).
#' @param pv.threshold P-value threshold for IV selection (e.g., 5e-8).
#' @param kappa_thres Numeric threshold for condition number filtering. Default is 10.
#' @param warnings Logical. Print data standardization reminder? Default is TRUE.
#'
#' @export
RaoBlackwellCorrect <- function(BETA_Select, SE_Select, Rxy,
                                gcov = matrix(0, nrow = ncol(Rxy), ncol = ncol(Rxy)),
                                ldsc = rep(0, nrow(BETA_Select)),
                                eta = 0.5,
                                pv.threshold,
                                kappa_thres = 10,
                                warnings = TRUE) {

  if (warnings) {
    cat("Please standardize data such that BETA = Zscore/sqrt(n) and SE = 1/sqrt(n)\n")
  }

  # Input validation
  stopifnot(all(dim(BETA_Select) == dim(SE_Select)))
  stopifnot(ncol(BETA_Select) == nrow(Rxy))
  stopifnot(nrow(Rxy) == ncol(Rxy))

  BETA_Select <- as.matrix(BETA_Select)
  SE_Select <- as.matrix(SE_Select)
  Rxy <- as.matrix(Rxy)

  m  <- nrow(BETA_Select)
  p1 <- ncol(BETA_Select)
  p  <- p1 - 1
  cutoff <- qchisq(pv.threshold, p, lower.tail = FALSE)

  # Format names
  if (is.null(rownames(BETA_Select))) {
    rownames(BETA_Select) <- paste0("SNP", seq_len(m))
  }
  if (is.null(colnames(BETA_Select))) {
    colnames(BETA_Select) <- c(paste0("Exp", seq_len(p)), "Outcome")
  }
  rownames(SE_Select) <- rownames(BETA_Select)
  colnames(SE_Select) <- colnames(BETA_Select)

  # Initialize storage
  BETA_RB <- BETA_Select
  SE_RB   <- matrix(NA_real_, m, p1)
  COV_RB  <- vector("list", m)
  corrected_initial <- integer(0)

  # Batch processing
  for (i in seq_len(m)) {
    sei <- SE_Select[i, ]
    Sigma_full <- t(t(Rxy) * sei) * sei + ldsc[i] * gcov

    out <- rb_onlyexposure_analytic_one(
      beta_hat = BETA_Select[i, ],
      Sigma_full = Sigma_full,
      eta = eta,
      cutoff = cutoff,
      p = p
    )

    # Accept if P_accept is valid and >= 0.1
    if (is.finite(out$P_accept) && (out$P_accept >= 0.1)) {
      BETA_RB[i, ] <- out$beta_RB
      SE_RB[i, ]   <- out$se_RB
      COV_RB[[i]]  <- out$cov_RB
      corrected_initial <- c(corrected_initial, i)
    } else {
      # Fallback for unaccepted SNPs
      BETA_RB[i, ] <- BETA_Select[i, ]
      SE_RB[i, ]   <- NA_real_
      M <- matrix(0, p1, p1)
      M[p1, p1] <- -1
      COV_RB[[i]] <- M
    }
  }

  # Condition Number Filtering (kappa_thres)
  ind_final <- integer(0)
  for (orig_i in corrected_initial) {
    if (is.na(SE_RB[orig_i, p1]) || SE_RB[orig_i, p1] <= 0) next

    sevec <- SE_Select[orig_i, ]
    kapparxy <- kappa(t(t(Rxy) * sevec) * sevec + ldsc[orig_i] * gcov)
    S <- COV_RB[[orig_i]]

    if (all(diag(S) > 0) && kappa(S) <= (kappa_thres * kapparxy)) {
      ind_final <- c(ind_final, orig_i)
    }
  }

  if (length(ind_final) == 0) {
    warning("No SNPs passed correction and filtering.")
    return(list(
      bX_RB = matrix(numeric(0), 0, p),
      bXse_RB = matrix(numeric(0), 0, p),
      by_RB = numeric(0),
      byse_RB = numeric(0),
      COV_RB = list(),
      CORRECTED_INDICES0 = corrected_initial,
      CORRECTED_INDICES = integer(0)
    ))
  }

  # Extract final results
  bX_RB <- BETA_RB[ind_final, 1:p, drop = FALSE]
  bXse_RB <- SE_RB[ind_final, 1:p, drop = FALSE]
  by_RB <- BETA_RB[ind_final, p + 1]
  byse_RB <- SE_RB[ind_final, p + 1]
  COV_RB_final <- COV_RB[ind_final]

  # Formatting output names for p=1 case
  if (p == 1) {
    bX_RB <- as.vector(bX_RB)
    bXse_RB <- as.vector(bXse_RB)
    names(bX_RB) <- names(bXse_RB) <- rownames(BETA_Select)[ind_final]
  } else {
    rownames(bX_RB) <- rownames(bXse_RB) <- rownames(BETA_Select)[ind_final]
    colnames(bX_RB) <- colnames(bXse_RB) <- colnames(BETA_Select)[1:p]
  }
  names(by_RB) <- names(byse_RB) <- rownames(BETA_Select)[ind_final]

  return(list(
    bX_RB = bX_RB,
    bXse_RB = bXse_RB,
    by_RB = by_RB,
    byse_RB = byse_RB,
    COV_RB = COV_RB_final,
    CORRECTED_INDICES0 = corrected_initial,
    CORRECTED_INDICES = ind_final
  ))
}

# ============================================================
# Internal helper functions
# ============================================================

g1g2_from_pchisq <- function(df, kappa, t) {
  stopifnot(df >= 1, kappa >= 0, t >= 0)

  logP0 <- pchisq(t, df = df,     ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP2 <- pchisq(t, df = df + 2, ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP4 <- pchisq(t, df = df + 4, ncp = kappa, lower.tail = FALSE, log.p = TRUE)

  if (!is.finite(logP0)) {
    return(list(P = 0, logP = logP0, g1 = NA_real_, g2 = NA_real_,
                r2 = NA_real_, r4 = NA_real_))
  }

  r2 <- exp(logP2 - logP0)
  r4 <- exp(logP4 - logP0)

  g1 <- 0.5 * (r2 - 1)
  P2_over_P <- 0.25 * (r4 - 2 * r2 + 1)
  g2 <- P2_over_P - g1^2

  list(P = exp(logP0), logP = logP0, g1 = g1, g2 = g2, r2 = r2, r4 = r4)
}

rb_onlyexposure_analytic_one <- function(beta_hat, Sigma_full, eta, cutoff, p = length(beta_hat) - 1) {
  beta_hat <- as.numeric(beta_hat)
  p1 <- length(beta_hat)
  stopifnot(is.matrix(Sigma_full), nrow(Sigma_full) == p1, ncol(Sigma_full) == p1)
  stopifnot(p >= 1, p < p1, eta > 0, cutoff >= 0)

  bx  <- beta_hat[1:p]
  Sxx <- Sigma_full[1:p, 1:p, drop = FALSE]

  cholS <- tryCatch(chol(Sxx), error = function(e) return(NULL))
  if (is.null(cholS)) {
    return(list(P_accept = 0)) # Failsafe instead of crash
  }
  Sxx_inv <- chol2inv(cholS)

  q     <- as.numeric(t(bx) %*% Sxx_inv %*% bx)
  kappa <- q / (eta^2)
  tval  <- cutoff / (eta^2)

  gg <- g1g2_from_pchisq(df = p, kappa = kappa, t = tval)
  Pacc <- gg$P
  if (!is.finite(Pacc) || Pacc <= 0) {
    return(list(P_accept = 0)) # Failsafe instead of 'stop()' crash
  }
  g1 <- gg$g1
  g2 <- gg$g2

  Ex   <- 2 * g1 * bx
  Covx <- eta^2 * (1 + 2 * g1) * Sxx + 4 * g2 * tcrossprod(bx)

  Sxy <- Sigma_full[1:p, (p+1):p1, drop = FALSE]
  Syx <- Sigma_full[(p+1):p1, 1:p, drop = FALSE]
  Syy <- Sigma_full[(p+1):p1, (p+1):p1, drop = FALSE]

  Amap <- Syx %*% Sxx_inv
  Er   <- as.numeric(Amap %*% Ex)

  Schur <- Syy - Amap %*% Sxy
  Covr  <- Amap %*% Covx %*% t(Amap) + eta^2 * Schur
  Covxr <- Covx %*% Sxx_inv %*% Sxy

  E_e <- c(Ex, Er)

  Cov_e <- matrix(0, p1, p1)
  Cov_e[1:p, 1:p] <- Covx
  Cov_e[(p+1):p1, (p+1):p1] <- Covr
  Cov_e[1:p, (p+1):p1] <- Covxr
  Cov_e[(p+1):p1, 1:p] <- t(Covxr)

  beta_RB <- beta_hat - E_e / (eta^2)
  cov_init_cond <- Cov_e / (eta^4)
  cov_RB <- (1 + 1/eta^2) * Sigma_full - cov_init_cond

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
