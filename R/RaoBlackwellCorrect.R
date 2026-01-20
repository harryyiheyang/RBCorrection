#' Rao-Blackwell Correction for Selected IVs in MVMR
#'
#' Performs Rao-Blackwell correction via rejection sampling for each selected IV
#' based on randomized multivariate chi-square test.
#'
#' @param BETA_Select A numeric matrix of dimension m x p: effect estimates for selected IVs. Outcome should be the last column.
#' @param SE_Select A numeric matrix of same dimension: standard errors. Outcome should be the last column.
#' @param Rxy A (p+1) x (p+1) numeric matrix: covariance matrix of Z statistics. Outcome should be the last column and row.
#' @param gcov A matrix (2 x 2) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.
#' @param eta Standard deviation of noise added to Z (default = 1).
#' @param pv.threshold P-value threshold used for selection.
#' @param B Number of samples used in rejection sampling (default = 1000).
#' @param onlyexposure A indicator of whether considering outcome when selecting IVs. Defaults to \code{TRUE}.
#' @param kappa_thres Numeric threshold controlling the removal of SNPs with unstable RBC covariance; if the condition number of the corrected covariance exceeds kappa_thres * kappa(Rxy) the SNP is discarded. Defaults to \code{10}.
#' @param warnings A indicator of whether printing warnings. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{BETA_RB}}{Corrected effect estimates (m x p).}
#'   \item{\code{SE_RB}}{Corrected standard errors (m x p).}
#' }
#'
#' @export
RaoBlackwellCorrect <- function(BETA_Select, SE_Select, Rxy, gcov=0*diag(BETA_Select[1,]), ldsc=0*BETA_Select[,1], eta = 0.5, pv.threshold,
                                B = 1000, kappa_thres = 10, onlyexposure = TRUE, warnings = TRUE, n_threads = min(1,parallel::detectCores()-2)) {
  if(warnings) {
    cat("Please standardize data such that BETA = Zscore/sqrt n and SE = 1/sqrt n\n")
  }

  stopifnot(all(dim(BETA_Select) == dim(SE_Select)))
  stopifnot(ncol(BETA_Select) == nrow(Rxy))
  stopifnot(nrow(Rxy) == ncol(Rxy))

  # Ensure matrix format with dimnames
  p = ncol(BETA_Select) - 1
  m = nrow(BETA_Select)

  if(onlyexposure == TRUE) {
    cutoff = qchisq(pv.threshold, p, lower.tail = FALSE)
  } else {
    cutoff = qchisq(pv.threshold, p + 1, lower.tail = FALSE)
  }

  BETA_Select <- as.matrix(BETA_Select)
  SE_Select <- as.matrix(SE_Select)
  Rxy <- as.matrix(Rxy)
  RxyList=RxysqrtList=array(0,c(p+1,p+1,m))
  for(i in 1:m){
  sei=SE_Select[i,]
  Rxyi=t(t(Rxy)*sei)*sei+ldsc[i]*gcov
  RxyList[,,i]=Rxyi
  RxysqrtList[,,i]=matrixsqrt(Rxyi)$w
  }

  if (is.null(rownames(BETA_Select))) rownames(BETA_Select) <- paste0("SNP", seq_len(nrow(BETA_Select)))
  if (is.null(colnames(BETA_Select))) colnames(BETA_Select) <- paste0("Exp", seq_len(ncol(BETA_Select)))
  rownames(SE_Select) <- rownames(BETA_Select)
  colnames(SE_Select) <- colnames(BETA_Select)

  res <- RaoBlackwell(
    beta_select = BETA_Select,
    se_select = SE_Select,
    RxyList = RxyList,
    RxysqrtList = RxysqrtList,
    eta = eta,
    cutoff = cutoff,
    B = B,
    min_accept = floor(B/2),
    onlyexposure = onlyexposure,
    n_threads = n_threads
  )

  for(i in res$CORRECTED_INDICES) {
    res$COV_RB[[i]] = RxyList[,,i]*(1 + 1/eta^2) - res$COV_RB[[i]]
    s = sqrt(diag(res$COV_RB[[i]]))
    s[is.na(s)] = 0
    res$SE_RB[i, ] = s
  }

  condition_check = function(res, thres = 3) {
    ind1 = which(!is.na(res$SE_RB[, p + 1]) & res$SE_RB[, p + 1] > 0)
    ind2 = integer(0)

    for (i in seq_along(res$COV_RB)) {
      sevec   = SE_Select[i, ]
      kapparxy = kappa(t(t(Rxy) * sevec) * sevec)+ldsc[i]*gcov
      S        = res$COV_RB[[i]]
      i1       = (diag(S) <= 0)
      i2_flag  = kappa(S) > (thres * kapparxy)

      if (sum(i1) + sum(i2_flag) == 0) {
        ind2 = c(ind2, i)
      }
    }

    ind = intersect(ind1, ind2)
    return(ind)
  }

  ind = condition_check(res, thres = kappa_thres)
  ind = intersect(res$CORRECTED_INDICES,ind)

  res$BETA_RB = res$BETA_RB[ind, ]
  res$SE_RB = res$SE_RB[ind, ]
  res$COV_RB = res$COV_RB[ind]

  rownames(res$BETA_RB) <- rownames(BETA_Select[ind, ])
  colnames(res$BETA_RB) <- colnames(BETA_Select)
  rownames(res$SE_RB) <- rownames(BETA_Select[ind, ])
  colnames(res$SE_RB) <- colnames(BETA_Select)

  bX_RB = res$BETA_RB[, 1:p]
  bXse_RB = res$SE_RB[, 1:p]
  by_RB = res$BETA_RB[, 1 + p]
  byse_RB = res$SE_RB[, 1 + p]

  if(p == 1) {
    bX_RB = as.vector(bX_RB)
    bXse_RB = as.vector(bXse_RB)
    names(bX_RB)=names(bXse_RB)=rownames(BETA_Select[ind, ])
  }

  return(list(bX_RB = bX_RB, bXse_RB = bXse_RB, by_RB = by_RB, byse_RB = byse_RB, COV_RB = res$COV_RB,
              CORRECTED_INDICES0 = res$CORRECTED_INDICES,
              CORRECTED_INDICES = ind))
}
