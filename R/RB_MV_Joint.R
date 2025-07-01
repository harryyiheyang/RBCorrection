#' Joint IV Selection and Rao-Blackwell Correction for Multivariable MR
#'
#' Performs joint IV selection using multivariate chi-square across exposures, followed by
#' Rao-Blackwell correction of selected SNPs in a simulation context.
#'
#' @param BETAMatrix A numeric matrix of SNP-exposure effect estimates (rows: SNPs, cols: exposures). Outcome should be the last column.
#' @param SEMatrix A numeric matrix of corresponding standard errors. Must match BETAMatrix in size and dimnames.Outcome should be the last column.
#' @param Rxy A (p+1 x p+1) covariance matrix of exposure effects. Outcome should be the last column and row.
#' @param P_threshold P-value threshold for selection (default = 0.05).
#' @param eta Standard deviation of noise for selection randomization (default = 1).
#' @param B Number of Monte Carlo samples for Rao-Blackwell correction (default = 10000).
#'
#' @return A list with two matrices: BETA_Sub (estimates) and SE_Sub (total SEs).
#' @export
RB_MV_Joint <- function(BETAMatrix, SEMatrix, Rxy,
                                  P_threshold = 0.05, eta = 1, B = 1000) {
  stopifnot(all(dim(BETAMatrix) == dim(SEMatrix)))
  if (is.null(rownames(BETAMatrix))) rownames(BETAMatrix) <- paste0("V", 1:nrow(BETAMatrix))
  if (is.null(rownames(SEMatrix))) rownames(SEMatrix) <- rownames(BETAMatrix)
  if (is.null(colnames(BETAMatrix))) colnames(BETAMatrix) <- paste0("Exp", 1:ncol(BETAMatrix))
  if (is.null(colnames(SEMatrix))) colnames(SEMatrix) <- colnames(BETAMatrix)

  m <- nrow(BETAMatrix)
  p <- ncol(BETAMatrix)-1
  cutoff <- qchisq(P_threshold, df = p, lower.tail = FALSE)

  ZMatrix <- BETAMatrix / SEMatrix
  Z_quasi <- ZMatrix + MASS::mvrnorm(n = m, mu = rep(0, p+1), Sigma = Rxy) * eta
  Thetaxx <- solve(Rxy[1:p,1:p])

  chisq_stat <- apply(Z_quasi[,1:p], 1, function(z) sum(z * matrixVectorMultiply(Thetaxx, z)))
  IV_select <- which(chisq_stat > cutoff)

  if (length(IV_select) == 0) {
    return(list(
      BETA_Sub = matrix(NA, nrow = 0, ncol = p+1, dimnames = list(NULL, colnames(BETAMatrix))),
      SE_Sub = matrix(NA, nrow = 0, ncol = p+1, dimnames = list(NULL, colnames(BETAMatrix)))
    ))
  }

  BETA_Select <- BETAMatrix[IV_select, , drop = FALSE]
  SE_Select <- SEMatrix[IV_select, , drop = FALSE]

  Resampling <- RaoBlackwellCorrect_R(
    BETA_Select = BETA_Select,
    SE_Select = SE_Select,
    Rxy = Rxy,
    eta = eta,
    cutoff = cutoff,
    B = B,
    onlyexposure=T
  )


  return(Resampling)
}
