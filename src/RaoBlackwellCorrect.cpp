// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <random>  // for thread-local RNG
using namespace Rcpp;

// [[Rcpp::export]]
List RaoBlackwell(NumericMatrix beta_select,
                         NumericMatrix se_select,
                         NumericMatrix Rxy,
                         NumericMatrix Rxysqrt,
                         double eta,
                         double cutoff,
                         int B,
                         bool onlyexposure = true,
                         int n_threads = 1) {
  int m = beta_select.nrow(), p1 = beta_select.ncol();
  int p = p1 - 1;

  arma::mat Beta(beta_select.begin(), m, p1, false);
  arma::mat SE(se_select.begin(), m, p1, false);
  arma::mat R(Rxy.begin(), p1, p1, false);

  if (!R.is_sympd()) stop("Rxy is not symmetric positive definite");

  arma::mat Rinv = arma::inv_sympd(R);
  arma::mat Rinv_sub;
  if (onlyexposure) {
    Rinv_sub = arma::inv_sympd(R.submat(0, 0, p - 1, p - 1));
  }

  arma::mat Rxysqrt_mat(Rxysqrt.begin(), Rxysqrt.nrow(), Rxysqrt.ncol(), false);
  arma::mat L = eta * Rxysqrt_mat;

  arma::mat BETA_RB(m, p1, arma::fill::zeros);
  arma::mat SE_RB(m, p1, arma::fill::zeros);
  std::vector<arma::mat> CovList(m);  // to store covariance matrices
  std::vector<int> corrected_indices;  // to store indices of corrected IVs

#ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < m; ++i) {
    // Thread-local random generator
    std::random_device rd;
    std::mt19937 rng(rd() + i * 9973 + omp_get_thread_num() * 99991);
    std::normal_distribution<double> norm(0.0, 1.0);

    arma::rowvec beta_i = Beta.row(i);
    arma::rowvec se_i = SE.row(i);

    std::vector<arma::rowvec> accepted;
    int accept_count = 0;

    for (int b = 0; b < 2 * B; ++b) {
      if (accept_count >= B) break;

      arma::vec z(p1);
      for (int j = 0; j < p1; ++j) {
        z(j) = norm(rng);
      }

      arma::vec e = L.t() * z;
      arma::vec z0 = arma::conv_to<arma::vec>::from(beta_i / se_i);
      arma::vec z_star = z0 + e;

      double chisq;
      if (onlyexposure) {
        arma::vec z_sub = z_star.subvec(0, p - 1);
        chisq = arma::as_scalar(z_sub.t() * Rinv_sub * z_sub);
      } else {
        chisq = arma::as_scalar(z_star.t() * Rinv * z_star);
      }

      if (chisq > cutoff) {
        arma::vec e_se = e % arma::trans(se_i);
        arma::rowvec beta_star = beta_i - arma::trans(e_se) / (eta * eta);
        accepted.push_back(beta_star);
        ++accept_count;
      }
    }

    if (accepted.size() >= 10) {
      // Sufficient samples: use correction
      arma::mat acc_mat(accepted.size(), p1);
      for (size_t j = 0; j < accepted.size(); ++j) {
        acc_mat.row(j) = accepted[j];
      }
      BETA_RB.row(i) = arma::mean(acc_mat, 0);
      SE_RB.row(i) = arma::stddev(acc_mat, 0, 0);
      CovList[i] = arma::cov(acc_mat, 0);

      // Record this IV as corrected (thread-safe)
#pragma omp critical
{
  corrected_indices.push_back(i + 1);  // R uses 1-based indexing
}

    } else {
      // Insufficient samples: use original data
      BETA_RB.row(i) = beta_i;
      SE_RB.row(i) = se_i;
      CovList[i] = arma::mat(p1, p1, arma::fill::zeros);
    }
  }

  NumericMatrix beta_out = wrap(BETA_RB);
  NumericMatrix se_out = wrap(SE_RB);
  beta_out.attr("dimnames") = beta_select.attr("dimnames");
  se_out.attr("dimnames") = se_select.attr("dimnames");

  List cov_list(CovList.size());
  for (int i = 0; i < m; ++i) {
    cov_list[i] = wrap(CovList[i]);
  }

  return List::create(Named("BETA_RB") = beta_out,
                      Named("SE_RB") = se_out,
                      Named("COV_RB") = cov_list,
                      Named("CORRECTED_INDICES") = corrected_indices);
}
