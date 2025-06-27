#' Rao-Blackwell Correction for Selected IVs in MVMR
#'
#' Performs Rao-Blackwell correction via rejection sampling for each selected IV
#' based on randomized multivariate chi-square test.
#'
#' @param BETA_Select A numeric matrix of dimension m x p: effect estimates for selected IVs. Outcome should be the last column.
#' @param SE_Select A numeric matrix of same dimension: standard errors. Outcome should be the last column.
#' @param Rxy A (p+1) x (p+1) numeric matrix: covariance matrix of Z statistics. Outcome should be the last column and row.
#' @param eta Standard deviation of noise added to Z (default = 1).
#' @param cutoff Chi-square threshold used for selection.
#' @param B Number of samples used in rejection sampling (default = 1000).
#' @param onlyexposure A indicator of whether considering outcome when selecting IVs. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{BETA_RB}}{Corrected effect estimates (m x p).}
#'   \item{\code{SE_RB}}{Corrected standard errors (m x p).}
#' }
#'
#' @export
RaoBlackwellCorrect_R <- function(BETA_Select, SE_Select, Rxy, eta = 1, cutoff,
                                  B = 1000, onlyexposure=T) {
  stopifnot(all(dim(BETA_Select) == dim(SE_Select)))
  stopifnot(ncol(BETA_Select) == nrow(Rxy))
  stopifnot(nrow(Rxy) == ncol(Rxy))

  # Ensure matrix format with dimnames
  BETA_Select <- as.matrix(BETA_Select)
  SE_Select <- as.matrix(SE_Select)
  Rxy <- as.matrix(Rxy)  # fix here

  # Ensure row and column names are present
  if (is.null(rownames(BETA_Select))) rownames(BETA_Select) <- paste0("SNP", seq_len(nrow(BETA_Select)))
  if (is.null(colnames(BETA_Select))) colnames(BETA_Select) <- paste0("Exp", seq_len(ncol(BETA_Select)))
  rownames(SE_Select) <- rownames(BETA_Select)
  colnames(SE_Select) <- colnames(BETA_Select)

  res <- .Call(`_RBCorrection_RaoBlackwellCorrect`, beta_select=BETA_Select, se_select=SE_Select,
               Rxy=Rxy, eta=eta,cutoff=cutoff, B=B, onlyexposure=onlyexposure,n_threads=floor(parallel::detectCores() / 2))
  res$SE_RB=sqrt(SE_Select^2*(1+1/eta^2))

  rownames(res$BETA_RB) <- rownames(BETA_Select)
  colnames(res$BETA_RB) <- colnames(BETA_Select)
  rownames(res$SE_RB) <- rownames(BETA_Select)
  colnames(res$SE_RB) <- colnames(BETA_Select)

  return(res)
}
