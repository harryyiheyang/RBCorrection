#' Rao-Blackwell Correction for Selected IVs in MVMR
#'
#' Performs Rao-Blackwell correction for instrumental variables using a closed-form
#' solution based on noncentral chi-square distribution theory.
#'
#' @param BETA_Select A numeric matrix of dimension m x (p+1): effect estimates.
#'   Columns 1:p are exposures, column p+1 is the outcome.
#' @param SE_Select A numeric matrix of same dimension: standard errors.
#' @param Rxy A (p+1) x (p+1) correlation matrix of Z-statistics from LDSC intercepts.
#' @param gcov A (p+1) x (p+1) matrix of per-SNP genetic covariances (heritability/M from LDSC).
#'   Default is zero matrix.
#' @param ldsc A numeric vector of length m: LD scores for each SNP. Default is zero vector.
#' @param eta Randomization parameter controlling noise level (default = 0.5).
#' @param pv.threshold P-value threshold for IV selection (e.g., 5e-8).
#' @param prob_thres Minimum selection probability threshold (default = 0.1).
#'   Variants with selection probability < prob_thres are excluded as low-quality IVs.
#' @param warnings Logical. Print data standardization reminder? Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{bX_RB}}{Corrected exposure effect estimates.}
#'   \item{\code{bXse_RB}}{Corrected exposure standard errors.}
#'   \item{\code{by_RB}}{Corrected outcome effect estimates.}
#'   \item{\code{byse_RB}}{Corrected outcome standard errors.}
#'   \item{\code{COV_RB}}{List of corrected covariance matrices.}
#'   \item{\code{CORRECTED_INDICES}}{Indices of SNPs passing correction.}
#' }
#'
#' @export
RaoBlackwellCorrect <- function(BETA_Select, SE_Select, Rxy,
                                gcov = matrix(0, nrow = ncol(Rxy), ncol = ncol(Rxy)),
                                ldsc = rep(0, nrow(BETA_Select)),
                                eta = 0.5,
                                pv.threshold,
                                prob_thres = 0.1,
                                warnings = TRUE) {

  if (warnings) {
    cat("Please standardize data such that BETA = Zscore/sqrt(n) and SE = 1/sqrt(n)\n")
  }

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

  if (is.null(rownames(BETA_Select))) {
    rownames(BETA_Select) <- paste0("SNP", seq_len(m))
  }
  if (is.null(colnames(BETA_Select))) {
    colnames(BETA_Select) <- c(paste0("Exp", seq_len(p)), "Outcome")
  }
  rownames(SE_Select) <- rownames(BETA_Select)
  colnames(SE_Select) <- colnames(BETA_Select)

  BETA_RB <- BETA_Select
  SE_RB   <- matrix(NA_real_, m, p1)
  COV_RB  <- vector("list", m)
  corrected <- integer(0)

  for (i in seq_len(m)) {
    sei <- SE_Select[i, ]
    Sigma_full <- t(t(Rxy) * sei) * sei + ldsc[i] * gcov

    out <- rb_onlyexposure_one(
      beta_hat = BETA_Select[i, ],
      Sigma_full = Sigma_full,
      eta = eta,
      cutoff = cutoff,
      p = p
    )

    if (is.finite(out$P_accept) && out$P_accept >= prob_thres) {
      BETA_RB[i, ] <- out$beta_RB
      SE_RB[i, ]   <- out$se_RB
      COV_RB[[i]]  <- out$cov_RB
      corrected <- c(corrected, i)
    } else {
      BETA_RB[i, ] <- BETA_Select[i, ]
      SE_RB[i, ]   <- NA_real_
      M <- matrix(0, p1, p1)
      M[p1, p1] <- -1
      COV_RB[[i]] <- M
    }
  }

  if (length(corrected) == 0) {
    warning("No SNPs passed correction. Consider relaxing prob_thres or pv.threshold.")
    return(list(
      bX_RB = matrix(numeric(0), 0, p),
      bXse_RB = matrix(numeric(0), 0, p),
      by_RB = numeric(0),
      byse_RB = numeric(0),
      COV_RB = list(),
      CORRECTED_INDICES = integer(0)
    ))
  }

  bX_RB <- BETA_RB[corrected, 1:p, drop = FALSE]
  bXse_RB <- SE_RB[corrected, 1:p, drop = FALSE]
  by_RB <- BETA_RB[corrected, p + 1]
  byse_RB <- SE_RB[corrected, p + 1]
  COV_RB_final <- COV_RB[corrected]

  if (p == 1) {
    bX_RB <- as.vector(bX_RB)
    bXse_RB <- as.vector(bXse_RB)
    names(bX_RB) <- names(bXse_RB) <- rownames(BETA_Select)[corrected]
  } else {
    rownames(bX_RB) <- rownames(bXse_RB) <- rownames(BETA_Select)[corrected]
    colnames(bX_RB) <- colnames(bXse_RB) <- colnames(BETA_Select)[1:p]
  }
  names(by_RB) <- names(byse_RB) <- rownames(BETA_Select)[corrected]

  return(list(
    bX_RB = bX_RB,
    bXse_RB = bXse_RB,
    by_RB = by_RB,
    byse_RB = byse_RB,
    COV_RB = COV_RB_final,
    CORRECTED_INDICES = corrected
  ))
}

# ============================================================
# Internal Functions
# ============================================================

g1g2_from_pchisq <- function(df, kappa, t) {
  stopifnot(df >= 1, kappa >= 0, t >= 0)

  logP0 <- pchisq(t, df = df,     ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP2 <- pchisq(t, df = df + 2, ncp = kappa, lower.tail = FALSE, log.p = TRUE)
  logP4 <- pchisq(t, df = df + 4, ncp = kappa, lower.tail = FALSE, log.p = TRUE)

  if (!is.finite(logP0)) {
    return(list(P = 0, g1 = NA_real_, g2 = NA_real_))
  }

  r2 <- exp(logP2 - logP0)
  r4 <- exp(logP4 - logP0)

  g1 <- 0.5 * (r2 - 1)
  g2 <- 0.25 * (r4 - 2 * r2 + 1) - g1^2

  list(P = exp(logP0), g1 = g1, g2 = g2)
}

rb_onlyexposure_one <- function(beta_hat, Sigma_full, eta, cutoff, p = length(beta_hat) - 1) {
  beta_hat <- as.numeric(beta_hat)
  p1 <- length(beta_hat)
  stopifnot(is.matrix(Sigma_full), nrow(Sigma_full) == p1, ncol(Sigma_full) == p1)
  stopifnot(p >= 1, p < p1, eta > 0, cutoff >= 0)

  bx  <- beta_hat[1:p]
  Sxx <- Sigma_full[1:p, 1:p, drop = FALSE]

  cholS <- tryCatch(chol(Sxx), error = function(e) return(NULL))
  if (is.null(cholS)) {
    return(list(P_accept = 0))
  }
  Sxx_inv <- chol2inv(cholS)

  q     <- as.numeric(t(bx) %*% Sxx_inv %*% bx)
  kappa <- q / (eta^2)
  tval  <- cutoff / (eta^2)

  gg <- g1g2_from_pchisq(df = p, kappa = kappa, t = tval)
  Pacc <- gg$P
  if (!is.finite(Pacc) || Pacc <= 0) {
    return(list(P_accept = 0))
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
  cov_RB <- (1 + 1/eta^2) * Sigma_full - Cov_e / (eta^4)
  cov_RB <- 0.5 * (cov_RB + t(cov_RB))

  list(
    beta_RB  = beta_RB,
    cov_RB   = cov_RB,
    se_RB    = sqrt(pmax(diag(cov_RB), 0)),
    P_accept = Pacc
  )
}
