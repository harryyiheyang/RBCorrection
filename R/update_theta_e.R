#' Rao-Blackwellized theta estimation
#'
#' Computes Rao-Blackwellized estimates of theta using a precomputed covariance array
#' and applies residual correction. This function wraps the C++ backend.
#'
#' @param by A numeric vector of outcome estimates.
#' @param bx A numeric vector of exposure estimates.
#' @param byse A numeric vector of standard errors for `by`.
#' @param ThetaList A 3D array of dimension m × 2 × 2 containing covariance matrices for each observation.
#' @param vartheta A numeric vector of length 2 specifying the projection direction (e.g., c(1, 0)).
#' @param indvalid An integer vector (1-based R indices) specifying valid entries to be used for estimating theta.
#' @param n_threads Number of threads for parallel computation. If NULL, uses half of available cores.
#'
#' @return A list with elements:
#' \describe{
#'   \item{theta}{Estimated theta based on the valid subset.}
#'   \item{e}{A numeric vector of absolute residuals: \code{abs(by - hatx * theta) / byse}.}
#'   \item{loss}{Sum of squared residuals for the baseline model: \code{sum((abs(by - bx * theta) / byse)^2)}.}
#' }
#' @export
update_theta_e <- function(by, bx, byse, ThetaList, vartheta, indvalid, n_threads = NULL) {
  stopifnot(length(by) == length(bx),
            length(byse) == length(by),
            is.array(ThetaList),
            dim(ThetaList)[2] == 2,
            dim(ThetaList)[3] == 2,
            length(vartheta) == 2)

  m <- length(by)
  if (!identical(dim(ThetaList)[1], m)) stop("ThetaList must have dim[1] == length(by)")

  # Flatten 3D ThetaList into row-major vector as expected by C++ backend
  Theta_flat <- c(ThetaList[,1,1], ThetaList[,1,2], ThetaList[,2,1], ThetaList[,2,2])

  if (is.null(n_threads)) {
    n_threads <- max(1L, floor(parallel::detectCores() / 2))
  }

  .Call(`_RBCorrection_update_theta_e`, by, bx, byse, Theta_flat, vartheta, indvalid, n_threads)
}
