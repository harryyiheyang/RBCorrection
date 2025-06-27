// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::export]]
List update_theta_e(NumericVector by,
                                         NumericVector bx,
                                         NumericVector byse,
                                         NumericVector ThetaList,
                                         NumericVector vartheta,
                                         IntegerVector indvalid,
                                         int n_threads = 0) {

  int m = by.size();

  // Convert to Armadillo vectors (no copy)
  arma::vec bY(by.begin(), m, false);
  arma::vec bX(bx.begin(), m, false);
  arma::vec bYse(byse.begin(), m, false);
  arma::vec vtheta(vartheta.begin(), 2, false);

  // Extract theta matrices more efficiently
  arma::vec Omega00(ThetaList.begin(), m, false);           // First m elements
  arma::vec Omega01(ThetaList.begin() + m, m, false);       // Next m elements
  arma::vec Omega10(ThetaList.begin() + 2*m, m, false);     // Next m elements
  arma::vec Omega11(ThetaList.begin() + 3*m, m, false);     // Last m elements

  const double vtheta0 = vtheta[0];
  const double vtheta1 = vtheta[1];

  // Vectorized computation of h for all observations
  arma::vec h = vtheta0 * (Omega00 * vtheta0 + Omega10 * vtheta1) +
    vtheta1 * (Omega01 * vtheta0 + Omega11 * vtheta1);

  // Vectorized computation of numerator
  arma::vec numerator = bX % (Omega00 * vtheta0 + Omega01 * vtheta1) +
    bY % (Omega10 * vtheta0 + Omega11 * vtheta1);

  // Vectorized tildex computation
  arma::vec hatx = numerator / h;

  // Vectorized g and H computation
  arma::vec g = hatx % ((bX - hatx) % Omega01 + bY % Omega11);
arma::vec H = arma::square(hatx) % Omega11;

// Convert indvalid to 0-based indexing
arma::uvec valid_idx(indvalid.size());
for (int i = 0; i < indvalid.size(); ++i) {
  valid_idx[i] = indvalid[i] - 1;
}

// Vectorized theta computation
double theta = arma::sum(g(valid_idx)) / arma::sum(H(valid_idx));

// Fully vectorized computation of e and e1 (like R's abs() function)
  arma::vec e = arma::abs(bY - hatx * theta) / bYse;
  arma::vec e1 = arma::abs(bY - bX * theta) / bYse;

  // Vectorized loss computation (like R's sum(e1[valid]^2))
double loss = arma::sum(arma::square(e1(valid_idx)));

return List::create(
  Named("theta") = theta,
  Named("e") = as<NumericVector>(wrap(e)),
  Named("loss") = loss
);
}
