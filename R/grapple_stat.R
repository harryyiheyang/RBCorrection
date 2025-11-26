grapple_stat <- function(RxyList, theta, e, n_threads = 1) {
p <- length(theta)
if (p == 1L) {
out <- grapple_stat_uni_cpp(RxyList, as.numeric(theta), e, n_threads = n_threads)
} else {
out <- grapple_stat_multi_cpp(RxyList, as.numeric(theta), e, n_threads = n_threads)
out$var_vec=as.vector(out$var_vec)
}
out
}
