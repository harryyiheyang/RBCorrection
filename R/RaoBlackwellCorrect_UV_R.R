#' Rao-Blackwell Correction via Resampling for Univariate MR
#'
#' Performs Rao-Blackwellized correction of SNP-exposure estimates via rejection sampling
#' using randomized Z-scores. This version wraps a fast C++ implementation.
#'
#' @param gamma A numeric vector of effect estimates (e.g. SNP-exposure).
#' @param sigma A numeric vector of standard errors, same length as \code{gamma}.
#' @param cutoff A numeric threshold on the randomized Z-statistic (i.e., abs(Z + noise) > cutoff).
#' @param B Number of Monte Carlo samples to use in rejection sampling (default = 1000).
#' @param eta Standard deviation of added noise (default = 1).
#' @param n_threads Number of threads to use (default = parallel::detectCores() / 2).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{\code{BETA_RB}}{Rao-Blackwell corrected effect estimates.}
#'   \item{\code{SE_RB}}{Estimated standard errors based on resampling.}
#' }
#'
#' @export
RaoBlackwellCorrect_UV_R <- function(gamma, sigma, cutoff, B = 1000, eta = 1,
                                     n_threads = round(parallel::detectCores() / 2)) {
  stopifnot(length(gamma) == length(sigma))
  stopifnot(is.numeric(gamma), is.numeric(sigma))

  res <- .Call(`_RBCorrection_RaoBlackwellCorrect_UV`, as.numeric(gamma), as.numeric(sigma),
               as.numeric(cutoff), as.integer(B), as.numeric(eta), as.integer(n_threads))

  out <- data.frame(BETA_RB = res$BETA_RB, SE_RB = sqrt(res$SE_RB^2+sigma^2))
  return(out)
}
