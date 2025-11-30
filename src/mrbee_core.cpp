// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//' Sum selected slices of a 3D array (cube) with OpenMP
 //'
 //' @param RxyList A (p x p x m) or (p+1 x p+1 x m) array, treated as arma::cube.
 //' @param indvalid Integer indices (1-based, as in R) indicating which slices
 //'   to sum over.
 //' @param n_threads Number of OpenMP threads to use (default 1).
 //'
 //' @return A p x p (or (p+1 x p+1)) matrix equal to
 //'   \code{Reduce("+", RxyList[,,indvalid])}.
 //'
 // [[Rcpp::export]]
 arma::mat biasterm(const arma::cube &RxyList,
                    const arma::uvec &indvalid,
                    int n_threads = 1) {

   int n_rows   = RxyList.n_rows;
   int n_cols   = RxyList.n_cols;
   int n_slices = RxyList.n_slices;

   arma::mat X(n_rows, n_cols, arma::fill::zeros);

#ifdef _OPENMP
   if (n_threads > 0)
     omp_set_num_threads(n_threads);

   //
#pragma omp parallel
{
  arma::mat X_local(n_rows, n_cols, arma::fill::zeros);

#pragma omp for nowait
  for (unsigned int k = 0; k < indvalid.n_elem; ++k) {
    int idx = indvalid[k] - 1;  // R (1-based) -> C++ (0-based)
    if (idx >= 0 && idx < n_slices) {
      X_local += RxyList.slice(idx);
    } else {
    }
  }

#pragma omp critical
{
  X += X_local;
}
}
#else
for (unsigned int k = 0; k < indvalid.n_elem; ++k) {
  int idx = indvalid[k] - 1;
  if (idx >= 0 && idx < n_slices) {
    X += RxyList.slice(idx);
  } else {
    stop("Index in indvalid is out of bounds for RxyList.");
  }
}
#endif

return X;
 }
