#' Rao-Blackwellized Estimate using Analytic Formula
#' Reference: Ma, Wang, Wu (2020)
#'
#' @param gamma Numeric vector of estimated effects (e.g., BETA)
#' @param sigma Numeric vector of standard errors (e.g., SE)
#' @param cutoff Selection threshold on z-score
#' @param eta Standard deviation of the injected noise
#'
#' @return data.frame with BETA_RB and SE_RB
#' @export
RB_UV_Analytic <- function(gamma, sigma, cutoff, eta = 1.0) {
  stopifnot(length(gamma) == length(sigma))

  # Step 1: randomized selection
  z_base <- gamma / sigma
  e <- rnorm(length(gamma), mean = 0, sd = eta)
  z_rand <- z_base + e
  sel <- abs(z_rand) > cutoff

  if (sum(sel) == 0) {
    return(data.frame(BETA_RB = numeric(0), SE_RB = numeric(0), select = integer(0)))
  }

  # Step 2: filter variables
  z_sel <- z_base[sel]
  gamma_sel <- gamma[sel]
  sigma_sel <- sigma[sel]

  # Step 3: apply formula
  Aj_plus <- -z_sel / eta + cutoff / eta
  Aj_minus <- -z_sel / eta - cutoff / eta

  phi_plus <- dnorm(Aj_plus)
  phi_minus <- dnorm(Aj_minus)
  Phi_plus <- pnorm(Aj_plus)
  Phi_minus <- pnorm(Aj_minus)

  denom <- 1 - Phi_plus + Phi_minus
  denom <- pmax(denom, 1e-12)  # 防止除零

  adjustment <- sigma_sel / eta * (phi_plus - phi_minus) / denom
  gamma_rb <- gamma_sel - adjustment

  part1 <- 1 - (1 / eta^2) * (Aj_plus * phi_plus - Aj_minus * phi_minus) / denom
  part2 <- (1 / eta^2) * ((phi_plus - phi_minus)^2) / (denom^2)
  var_rb <- sigma_sel^2 * (part1 + part2)
  se_rb <- sqrt(pmax(var_rb, 0))

  return(data.frame(BETA_RB = gamma_rb, SE_RB = se_rb, IVselect = which(sel)))
}
