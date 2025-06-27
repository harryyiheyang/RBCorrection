#' Joint Rao-Blackwellized Univariate MR Estimation for Simulation
#'
#' Performs joint Rao-Blackwell correction on simulated SNP-exposure and SNP-outcome associations
#' based on randomized instrument selection and rejection sampling.
#'
#' @param GWAS_X A data.frame with columns: SNP, BETA, SE (exposure).
#' @param GWAS_Y A data.frame with columns: SNP, BETA, SE (outcome), must match GWAS_X.
#' @param Rxy A 2x2 covariance matrix for the standardized Z-scores of exposure and outcome.
#' @param P_threshold P-value threshold for IV selection (default = 0.05).
#' @param eta Standard deviation of noise added for randomization (default = 1).
#' @param B Number of Monte Carlo samples for Rao-Blackwellization (default = 10000).
#' @importFrom MASS mvrnorm
#' @importFrom stats qchisq rnorm dnorm pnorm
#' @return A list of data.frames: GWAS_X_RB and GWAS_Y_RB, each with SNP, BETA_RB, SE_RB.
#' @export
RB_UV_Simulation_Joint <- function(GWAS_X, GWAS_Y, Rxy,
                                   P_threshold = 0.05, eta = 1, B = 1000) {
  stopifnot(all(c("SNP", "BETA", "SE") %in% names(GWAS_X)))
  stopifnot(all(c("SNP", "BETA", "SE") %in% names(GWAS_Y)))
  stopifnot(nrow(GWAS_X) == nrow(GWAS_Y))
  stopifnot(all(GWAS_X$SNP == GWAS_Y$SNP))
  stopifnot(all(dim(Rxy) == c(2, 2)))

  cutoff <- sqrt(qchisq(P_threshold, df = 1, lower.tail = FALSE))
  n <- nrow(GWAS_X)

  # Standardize Z and add bivariate Gaussian noise
  z_base_X <- GWAS_X$BETA / GWAS_X$SE
  z_base_Y <- GWAS_Y$BETA / GWAS_Y$SE
  Z_base <- cbind(z_base_X, z_base_Y)

  noise <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Rxy * eta^2)
  Z_rand <- Z_base + noise

  # Selection based on randomized Z for exposure
  sel <- abs(Z_rand[, 1]) > cutoff
  if (sum(sel) == 0) {
    return(list(
      GWAS_X_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric()),
      GWAS_Y_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric())
    ))
  }

  gamma <- GWAS_X$BETA[sel]
  Gamma <- GWAS_Y$BETA[sel]
  sigma_X <- GWAS_X$SE[sel]
  sigma_Y <- GWAS_Y$SE[sel]
  SNPs <- GWAS_X$SNP[sel]

  RB <- RaoBlackwellCorrect_R(
    BETA_Select = cbind(gamma, Gamma),
    SE_Select = cbind(sigma_X, sigma_Y),
    Rxy = Rxy,
    eta = eta,
    cutoff = cutoff^2,
    B = B,
    onlyexposure=T
  )

  out_X <- data.frame(SNP = SNPs, BETA_RB = RB$BETA_RB[, 1], SE_RB = RB$SE_RB[, 1])
  out_Y <- data.frame(SNP = SNPs, BETA_RB = RB$BETA_RB[, 2], SE_RB = RB$SE_RB[, 2])

  return(list(GWAS_X_RB = out_X, GWAS_Y_RB = out_Y, S_RB=RB$COV_RB))
}
