// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
List RaoBlackwell(NumericMatrix beta_select,
                  NumericMatrix se_select,
                  arma::cube RxyList,
                  arma::cube RxysqrtList,
                  double eta,
                  double cutoff,
                  int B,
                  int min_accept = 100,
                  bool onlyexposure = true,
                  int n_threads = 1) {

  int m = beta_select.nrow();
  int p1 = beta_select.ncol();
  int p = p1 - 1;

  arma::mat Beta(beta_select.begin(), m, p1, false);
  arma::mat SE(se_select.begin(), m, p1, false);

  arma::mat BETA_RB(m, p1, arma::fill::zeros);
  arma::mat SE_RB(m, p1, arma::fill::zeros);
  std::vector<arma::mat> CovList(m);
  std::vector<int> corrected_indices;

  const int max_draws = 5 * B;

#ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < m; ++i) {
    std::random_device rd;
    unsigned int seed = rd() + i * 9973;
#ifdef _OPENMP
    seed += omp_get_thread_num() * 99991;
#endif
    std::mt19937 rng(seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    arma::rowvec beta_i = Beta.row(i);

    arma::mat Sigma_i = RxyList.slice(i);
    arma::mat Sqrt_i  = RxysqrtList.slice(i);

    arma::mat Sigma_inv;
    if (onlyexposure) {
      Sigma_inv = arma::inv_sympd(Sigma_i.submat(0, 0, p - 1, p - 1));
    } else {
      Sigma_inv = arma::inv_sympd(Sigma_i);
    }

    arma::mat L = eta * Sqrt_i;

    arma::mat Z_raw(max_draws, p1);
    for (int draw = 0; draw < max_draws; ++draw) {
      for (int j = 0; j < p1; ++j) {
        Z_raw(draw, j) = norm(rng);
      }
    }
    Z_raw.each_row() -= arma::mean(Z_raw, 0);

    arma::mat E = Z_raw * L.t();

    std::vector<arma::rowvec> accepted;
    accepted.reserve(B);

    for (int draw = 0; draw < max_draws; ++draw) {
      if ((int)accepted.size() >= B) break;

      arma::vec e = E.row(draw).t();
      arma::vec beta_aug = arma::trans(beta_i) + e;

      double chisq;
      if (onlyexposure) {
        arma::vec beta_aug_sub = beta_aug.subvec(0, p - 1);
        chisq = arma::as_scalar(beta_aug_sub.t() * Sigma_inv * beta_aug_sub);
      } else {
        chisq = arma::as_scalar(beta_aug.t() * Sigma_inv * beta_aug);
      }

      if (chisq > cutoff) {
        arma::rowvec beta_init = beta_i - arma::trans(e) / (eta * eta);
        accepted.push_back(beta_init);
      }
    }

    if ((int)accepted.size() >= min_accept) {
      arma::mat acc_mat(accepted.size(), p1);
      for (size_t j = 0; j < accepted.size(); ++j) {
        acc_mat.row(j) = accepted[j];
      }

      BETA_RB.row(i) = arma::mean(acc_mat, 0);
      SE_RB.row(i)   = arma::stddev(acc_mat, 0, 0);
      CovList[i]     = arma::cov(acc_mat, 0);

#pragma omp critical
{
  corrected_indices.push_back(i + 1);
}
    } else {
      BETA_RB.row(i) = beta_i;
      SE_RB.row(i)   = SE.row(i);

      CovList[i] = arma::mat(p1, p1, arma::fill::zeros);
      CovList[i](p1 - 1, p1 - 1) = -1.0;
    }
  }

  NumericMatrix beta_out = wrap(BETA_RB);
  NumericMatrix se_out = wrap(SE_RB);

  beta_out.attr("dimnames") = beta_select.attr("dimnames");
  se_out.attr("dimnames") = se_select.attr("dimnames");

  List cov_list(m);
  for (int i = 0; i < m; ++i) {
    cov_list[i] = wrap(CovList[i]);
  }

  return List::create(
    Named("BETA_RB") = beta_out,
    Named("SE_RB") = se_out,
    Named("COV_RB") = cov_list,
    Named("CORRECTED_INDICES") = corrected_indices
  );
}
