// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>  // for std::max

using namespace Rcpp;

//' GRAPPLE stat (multivariable) in C++ with OpenMP
 //'
 //' @param RxyList A (p+1) x (p+1) x m cube, the k-th slice is the
 //'   (p+1) x (p+1) matrix for SNP k.
 //' @param theta A length-p vector of current causal effect estimates.
 //' @param e A length-m vector of residuals (so that e[i]^2 is r_i^2
 //'   in the Hessian derivation).
 //' @param n_threads Number of OpenMP threads.
 //'
 //' @return List with var_vec (length m), var_cor (m x p),
 //'   bias_correction (p x p) using the empirical Hessian with e^2.
 //'
 //' @export
 // [[Rcpp::export]]
 Rcpp::List grapple_stat_multi_cpp(const arma::cube& RxyList,
                                   const arma::vec& theta,
                                   const arma::vec& e,
                                   int n_threads = 1) {

   int p    = theta.n_elem;       // number of exposures
   int m    = RxyList.n_slices;   // number of SNPs
   int dim1 = RxyList.n_rows;     // should be p+1
   int dim2 = RxyList.n_cols;     // should be p+1

   if (dim1 != p + 1 || dim2 != p + 1) {
     Rcpp::stop("RxyList must have dimension (p+1, p+1, m) in grapple_stat_multi_cpp.");
   }
   if (e.n_elem != (unsigned int)m) {
     Rcpp::stop("Length of e must equal number of slices m in RxyList.");
   }

   // vartheta = (theta, -1)
   arma::vec vartheta(p + 1);
   vartheta.head(p) = theta;
   vartheta(p)      = -1.0;

   arma::vec var_vec(m);        // length m
   arma::mat var_cor(m, p);     // m x p

   // thread-local bias_cor: p x p x n_threads
   int n_th = std::max(n_threads, 1);
   arma::cube bias_thread(p, p, n_th, arma::fill::zeros);

#ifdef _OPENMP
   if (n_threads > 0) {
     omp_set_num_threads(n_threads);
   } else {
     n_threads = 1;
   }
#else
   n_threads = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int i = 0; i < m; ++i) {
     int tid = 0;
#ifdef _OPENMP
     tid = omp_get_thread_num();
#endif

     arma::mat Mi = RxyList.slice(i);  // (p+1) x (p+1)

     // var_vec[i] = vartheta' * Mi * vartheta
     double v_i = arma::as_scalar(vartheta.t() * Mi * vartheta);
     var_vec(i) = v_i;

     // hi = Mi[1:p,1:p] * theta - Mi[1:p,p+1]
     arma::mat Mxx = Mi.submat(0, 0, p-1, p-1);   // p x p
     arma::vec Mxy = Mi.submat(0, p, p-1, p);     // p x 1
     arma::vec hi  = Mxx * theta - Mxy;

     var_cor.row(i) = hi.t();

     // ---- 严格 Hessian：把原来 E[r_i^2]=1 的近似，改成观测 r_i^2 = e(i)^2 ----
     double r2 = e(i) * e(i);         // r_i^2
     double w  = r2 / (v_i * v_i);    // 对应 1 / v_i^2 * r_i^2
     bias_thread.slice(tid) += Mxx * w;
   }

   // sum bias_cor over threads
   arma::mat bias_cor(p, p, arma::fill::zeros);
   for (int t = 0; t < n_th; ++t) {
     bias_cor += bias_thread.slice(t);
   }

   return Rcpp::List::create(
     Rcpp::Named("var_vec")         = var_vec,
     Rcpp::Named("var_cor")         = var_cor,
     Rcpp::Named("bias_correction") = bias_cor
   );
 }

 //' GRAPPLE stat (univariable) in C++ with OpenMP
 //'
 //' @param RxyList A 2 x 2 x m cube. The k-th slice is the 2 x 2
 //'   matrix for SNP k.
 //' @param theta Scalar causal effect estimate.
 //' @param e A length-m vector of residuals (so that e[i]^2 is r_i^2
 //'   in the Hessian derivation).
 //' @param n_threads Number of OpenMP threads.
 //'
 //' @return List with var_vec (length m), var_cor (length m),
 //'   bias_correction (scalar) using the empirical Hessian with e^2.
 //'
 //' @export
 // [[Rcpp::export]]
 Rcpp::List grapple_stat_uni_cpp(const arma::cube& RxyList,
                                 double theta,
                                 const arma::vec& e,
                                 int n_threads = 1) {

   int m    = RxyList.n_slices;   // number of SNPs
   int dim1 = RxyList.n_rows;     // should be 2
   int dim2 = RxyList.n_cols;     // should be 2

   if (dim1 != 2 || dim2 != 2) {
     Rcpp::stop("RxyList must have dimension (2, 2, m) in grapple_stat_uni_cpp.");
   }
   if (e.n_elem != (unsigned int)m) {
     Rcpp::stop("Length of e must equal number of slices m in RxyList.");
   }

   // vartheta = (theta, -1)  (length 2)
   arma::vec vartheta(2);
   vartheta(0) = theta;
   vartheta(1) = -1.0;

   arma::vec var_vec(m);    // length m
   arma::vec var_cor(m);    // length m

   // thread-local bias scalar
   int n_th = std::max(n_threads, 1);
   arma::vec bias_thread(n_th, arma::fill::zeros);

#ifdef _OPENMP
   if (n_threads > 0) {
     omp_set_num_threads(n_threads);
   } else {
     n_threads = 1;
   }
#else
   n_threads = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int i = 0; i < m; ++i) {
     int tid = 0;
#ifdef _OPENMP
     tid = omp_get_thread_num();
#endif

     arma::mat Mi = RxyList.slice(i);  // 2 x 2

     // var_vec[i] = vartheta' * Mi * vartheta
     double v_i = arma::as_scalar(vartheta.t() * Mi * vartheta);
     var_vec(i) = v_i;

     // hi = Mi[1,1]*theta - Mi[1,2]  (0-based index: (0,0) and (0,1))
     double hi = Mi(0,0) * theta - Mi(0,1);
     var_cor(i) = hi;

     // 严格 Hessian：Mi(1,1) * r_i^2 / v_i^2
     double r2 = e(i) * e(i);
     double w  = r2 / (v_i * v_i);
     bias_thread(tid) += Mi(0,0) * w;
   }

   // sum bias over threads
   double bias_cor = arma::sum(bias_thread);

   return Rcpp::List::create(
     Rcpp::Named("var_vec")         = var_vec,
     Rcpp::Named("var_cor")         = var_cor,
     Rcpp::Named("bias_correction") = bias_cor
   );
 }
