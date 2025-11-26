// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//' MRcML_bXest in C++ with OpenMP (multivariable)
//'
//' @param ThetaList (p+1) x (p+1) x m array, as arma::cube
//' @param bX m x p matrix of original exposure effects
//' @param e length-m vector of residuals
//' @param theta length-p vector of causal effects
//' @param n_threads number of OpenMP threads
//' @return m x p matrix of bXest
//'
// [[Rcpp::export]]
arma::mat MRcML_bXest(const arma::cube &ThetaList,
                          const arma::mat  &bX,
                          const arma::vec  &e,
                          const arma::vec  &theta,
                          int n_threads = 1) {

  int m = bX.n_rows;
  int p = bX.n_cols;
  int p1 = ThetaList.n_rows;  // should be p+1

  if (ThetaList.n_slices != m)
    stop("ThetaList: number of slices must equal nrow(bX).");
  if (p1 != p + 1)
    stop("ThetaList: first dimension must be p+1.");

  arma::mat bXest = bX;

#ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#pragma omp parallel for
#endif
  for (int i = 0; i < m; ++i) {

    arma::mat Theta = ThetaList.slice(i);   // (p+1) x (p+1)

    // G: p x (p+1) = [I_p, theta]
    arma::mat G(p, p1, arma::fill::zeros);
    G.submat(0, 0, p-1, p-1) = arma::eye<arma::mat>(p, p);
    G.col(p) = theta;

    arma::mat GTheta = G * Theta;          // p x (p+1)
    arma::mat GTG    = GTheta * G.t();     // p x p

    // r = (bX[i,], e[i])
    arma::vec r(p1);
    for (int j = 0; j < p; ++j) r(j) = bX(i, j);
    r(p) = e(i);

    // bXest[i,] = solve(GTG) %*% GTheta %*% r
    arma::vec sol = arma::solve(GTG, GTheta * r);
    for (int j = 0; j < p; ++j) bXest(i, j) = sol(j);
  }

  return bXest;
}

//' MRcML_UV_bxest in C++ with OpenMP (univariable)
//'
//' @param ThetaList 2 x 2 x m array, as arma::cube
//' @param bx length-m vector of exposure effects
//' @param e  length-m vector of residuals
//' @param theta scalar causal effect
//' @param n_threads number of OpenMP threads
//' @return length-m vector of bxest
//'
// [[Rcpp::export]]
arma::vec MRcML_UV_bxest(const arma::cube &ThetaList,
                             const arma::vec  &bx,
                             const arma::vec  &e,
                             double theta,
                             int n_threads = 1) {

  int m  = bx.n_elem;
  int p1 = ThetaList.n_rows;  // should be 2

  if (ThetaList.n_slices != m)
    stop("ThetaList: number of slices must equal length(bx).");
  if (p1 != 2)
    stop("ThetaList: first dimension must be 2 (p+1 for p=1).");

  arma::vec bxest = bx;

  #ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  #pragma omp parallel for
  #endif
  for (int i = 0; i < m; ++i) {

    arma::mat Theta = ThetaList.slice(i);  // 2 x 2

    // G = c(1, theta)
    arma::vec G(2);
    G(0) = 1.0;
    G(1) = theta;

    // GTG = t(G) %*% Theta %*% G  (scalar)
    double GTG = arma::as_scalar(G.t() * Theta * G);

    // r = c(bx[i], e[i])
    arma::vec r(2);
    r(0) = bx(i);
    r(1) = e(i);

    // numerator = t(G) %*% Theta %*% r
    double num = arma::as_scalar(G.t() * Theta * r);

    bxest(i) = num / GTG;
  }

  return bxest;
}
