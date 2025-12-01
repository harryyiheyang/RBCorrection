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
                  NumericMatrix Rxy,
                  NumericMatrix Rxysqrt,
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
    arma::rowvec se_i = SE.row(i);

    arma::mat Z(max_draws, p1);
    for (int draw = 0; draw < max_draws; ++draw) {
      for (int j = 0; j < p1; ++j) {
        Z(draw, j) = norm(rng);
      }
    }

    Z.each_row() -= arma::mean(Z, 0);


    arma::mat E = Z * L.t();

    arma::vec z0 = arma::conv_to<arma::vec>::from(beta_i / se_i);

    std::vector<arma::rowvec> accepted;
    accepted.reserve(B);

    for (int draw = 0; draw < max_draws; ++draw) {
      if ((int)accepted.size() >= B) break;

      arma::vec e = E.row(draw).t();
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
      }
    }

    if ((int)accepted.size() >= min_accept) {
      arma::mat acc_mat(accepted.size(), p1);
      for (size_t j = 0; j < accepted.size(); ++j) {
        acc_mat.row(j) = accepted[j];
      }

      BETA_RB.row(i) = arma::mean(acc_mat, 0);
      SE_RB.row(i) = arma::stddev(acc_mat, 0, 0);
      CovList[i] = arma::cov(acc_mat, 0);

#pragma omp critical
{
  corrected_indices.push_back(i + 1);
}
    } else {
      BETA_RB.row(i) = beta_i;
      SE_RB.row(i) = se_i;

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
