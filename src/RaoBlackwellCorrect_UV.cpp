#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
List RaoBlackwellCorrect_UV(NumericVector gamma,
                            NumericVector sigma,
                            double cutoff,
                            int B,
                            double eta = 1.0,
                            int n_threads = 1) {
  int m = gamma.size();
  NumericVector beta_rb(m);
  NumericVector se_rb(m);

#ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif

#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    double beta_i = gamma[i];
    double se_i = sigma[i];

    std::vector<double> accepted;
    accepted.reserve(B / 2);  // heuristic

    for (int b = 0; b < B; ++b) {
      double eps = R::rnorm(0.0, eta);
      double z_rand = beta_i / se_i + eps;
      if (std::abs(z_rand) > cutoff) {
        double gamma_init = beta_i - se_i * eps / (eta * eta);
        accepted.push_back(gamma_init);
      }
    }

    int N = accepted.size();
    if (N > 0) {
      double sum = 0.0, sumsq = 0.0;
      for (int j = 0; j < N; ++j) {
        sum += accepted[j];
        sumsq += accepted[j] * accepted[j];
      }
      double mean = sum / N;
      beta_rb[i] = mean;
      se_rb[i] = std::sqrt(std::max(0.0, sumsq / N - mean * mean));
    } else {
      beta_rb[i] = NA_REAL;
      se_rb[i] = NA_REAL;
    }
  }

  return List::create(
    Named("BETA_RB") = beta_rb,
    Named("SE_RB") = se_rb
  );
}
