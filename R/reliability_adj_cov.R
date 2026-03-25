#' Reliability adjustment for Rao-Blackwell correction covariance matrices
#'
#' For each IV, limits the proportion of estimation error relative to total
#' variance to prevent unstable estimation caused by extremely weak instruments.
#'
#' @param bX A matrix (n x p) of GWAS effect sizes of p exposures.
#' @param Cov_RB A list of n matrices of Rao-Blackwell correction terms,
#'   each of dimension (p+1) x (p+1).
#' @param threshold Minimum reliability threshold. Defaults to 0.5.
#'
#' @return A list of n adjusted covariance matrices.
#' @export
reliability_adj_cov = function(bX, Cov_RB, threshold = 0.5) {
  n = nrow(bX)
  p = ncol(bX)
  Glist = Cov_RB
  for (i in 1:n) {
    G = Cov_RB[[i]]
    total.var = bX[i, ]^2
    error.var = diag(G[1:p, 1:p])
    reliability = (total.var - error.var) / total.var
    r = rep(1, p)
    ind = which(reliability < threshold)
    if (length(ind) > 0) {
      r[ind] = total.var[ind] / error.var[ind] * (1 - threshold)
    }
    r = c(sqrt(r), 1)
    Glist[[i]] = t(t(G) * r) * r
  }
  return(Glist)
}
