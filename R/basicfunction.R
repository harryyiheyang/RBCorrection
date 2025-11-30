
IVweight=function(byse,bXse,Rxy){
  bZse=cbind(bXse,byse)
  p=dim(bZse)[2]
  n=dim(bZse)[1]
  RxyList=array(0,c(n,p,p))
  for(i in 1:n){
    s=bZse[i,]
    RxyList[i,,]=t(t(Rxy)*s)*s
  }
  return(RxyList)
}

matrixsqrt=function(A){
  fit=matrixEigen(A)
  d=c(fit$value)
  d1=d*0
  d1[d>0]=1/d[d>0]
  d=sqrt(d)
  d1=sqrt(d1)
  A=matrixMultiply(fit$vector,t(fit$vector)*d)
  B=matrixMultiply(fit$vector,t(fit$vector)*d1)
  C=list(w=A,wi=B,eigenfit=fit)
  return(C)
}

imrpdetect=function(x,theta,RxyList,indvalid,var.est="robust",FDR=T,adjust.method="Sidak"){
  p=length(theta)
  if(var.est=="robust"){
    varx=stats::mad(x[indvalid])^2
  }
  if(var.est=="variance"){varx=stats::var(x[indvalid])}
  if(var.est=="ordinal"){
    varx=x*0
    for(i in 1:length(x)){
      varx[i]=c(RxyList[p+1,p+1,i]+t(theta)%*%RxyList[1:p,1:p,i]%*%theta-2*sum(theta*RxyList[p+1,1:p,i]))
    }
  }
  pv=stats::pchisq(x^2/varx,1,lower.tail=F)
  if(FDR==T){
    pv=FDRestimation::p.fdr(pvalues=pv,adjust.method=adjust.method,ties.method="random")$fdrs
  }
  return(as.vector(pv))
}

validadj <- function(vector1, vector2, tau) {
  diff <- length(vector2) / length(vector1)
  if (diff < tau) {
    missing_indices <- setdiff(1:length(vector1), vector2)
    sorted_missing_indices <- missing_indices[order(vector1[missing_indices])]
    num_to_add <- ceiling(tau * length(vector1)) - length(vector2)
    vector2 <- c(vector2, sorted_missing_indices[1:num_to_add])
  }
  return(vector2)
}

get_top_indices <- function(x, k) {
  x <- as.numeric(x)
  if (k == 0) {
    return(which(!is.na(x)))
  }
  x_sorted <- sort(x, decreasing = TRUE, na.last = NA)
  if (k >= length(x_sorted)) {
    return(integer(0))
  }
  threshold <- x_sorted[k + 1]
  which(x <= threshold & !is.na(x))
}

solve_tau <- function(tau_eq, tau_upper=5, n_grid = 30) {
  f0 <- tau_eq(0)
  fU <- tau_eq(tau_upper)
  if (f0 * fU < 0) {
    out <- uniroot(tau_eq, c(0, tau_upper),maxiter=100)
    return(out$root)
  }
  grid <- seq(0, tau_upper, length.out = n_grid + 1)
  vals <- vapply(grid, tau_eq, numeric(1))
  tau_hat <- grid[which.min(abs(vals))]
  return(tau_hat)
}

mcp = function(x, lam, a = 3) {
  b = abs(x)
  z = soft(x, lam) / (1 - 1 / a)
  z[which(b > a * lam)] = x[which(b > a * lam)]
  return(z)
}

pleio_adj <- function(gamma, max.prop.pleio = 0.5) {
n <- length(gamma)
nz_idx <- which(gamma != 0)
prop_nz <- length(nz_idx) / n
if (prop_nz <= max.prop.pleio) return(gamma)
abs_gamma <- abs(gamma[nz_idx])
thr <- quantile(abs_gamma, probs = 1 - max.prop.pleio)
gamma_new <- soft(gamma, thr/2)
return(gamma_new)
}


rho_mcp <- function(u, lambda = 2, gamma = 3) {
  u_abs <- abs(u)
  rho   <- u * 0
  ind1 <- u_abs <= lambda
  rho[ind1] <- 0.5 * u_abs[ind1]^2
  ind2 <- (u_abs > lambda) & (u_abs < gamma * lambda)
  if (any(ind2)) {
    rho[ind2] <- 0.5 * gamma * lambda^2 - 0.5 * (gamma * lambda - u_abs[ind2])^2 / (gamma - 1)
  }
  ind3 <- u_abs >= gamma * lambda
  if (any(ind3)) {
    rho[ind3] <- 0.5 * gamma * lambda^2
  }
  rho
}

rho_tukey <- function(r, c = 4.6851) {
  r_abs <- abs(r)
  rho   <- r * 0          # keep same shape as r
  ind1 <- r_abs <= c
  if (any(ind1)) {
    rc <- r[ind1] / c
    rho[ind1] <- (c^2 / 6) * (1 - (1 - rc^2)^3)
  }
  ind2 <- !ind1
  if (any(ind2)) {
    rho[ind2] <- c^2 / 6
  }
  rho
}
soft <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}
MCP_simulation=function(iter=10000,lambda=2,gamma=3){
  x=rnorm(iter)
  y=rho_mcp(x,lambda=lambda,gamma=gamma)
  return(mean(y))
}
