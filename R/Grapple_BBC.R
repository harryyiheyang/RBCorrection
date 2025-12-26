#' GRAPPLE with background bias correction (BBC).
#'
#' This function estimates causal effects using a robust GRAPPLE-type
#' likelihood method combined with background bias correction (BBC).
#' It uses an MCP-based loss to downweight outlying instruments and
#' accounts for sample overlap through LDSC and Rao-Blackwell corrections.
#'
#' @param by A vector (n x 1) of the GWAS effect sizes of the outcome
#'   after Rao-Blackwell correction.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures
#'   after Rao-Blackwell correction. If univariable, \code{bX} can be a vector
#'   of length n.
#' @param byse A vector (n x 1) of the standard errors of \code{by}
#'   after Rao-Blackwell correction.
#' @param bXse A matrix (n x p) of the standard errors of \code{bX}
#'   after Rao-Blackwell correction. If univariable, \code{bXse} can be a
#'   vector of length n.
#' @param cov_RB A list of length n, where each element is a (p+1 x p+1)
#'   matrix of Rao-Blackwell correction terms for the joint
#'   (exposures, outcome) GWAS estimates at a given SNP.
#' @param max.iter Maximum number of iterations for causal effect estimation.
#'   Defaults to 30.
#' @param max.eps Tolerance for the stopping criterion based on the change
#'   in the causal effect estimates. Defaults to 1e-4.
#' @param lambda Tuning parameter for the MCP loss used to downweight outliers.
#'   Defaults to 2.
#' @param a Concavity parameter of the MCP loss (often denoted gamma in the
#'   MCP literature). Must be greater than 1. Defaults to 3.
#' @param tau_upper Upper bound of the search range for the variance
#'   inflation parameter \code{tau} in the robust variance equation.
#'   Defaults to 10.
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   covariance (or variance) of the causal effect estimates. Defaults to 300.
#' @param n_threads The number of threads. Defaults to 4.
#' @param max.prop.pleio A numeric value in (0, 1). The maximum allowed proportion  of non-zero pleiotropic terms. Default is 0.5.
#' @param sampling.strategy "bootstrap" or "subsampling" (0.5 sample without replacement).

#' @return A list containing:
#' \describe{
#'   \item{theta}{A vector (length p) of the estimated causal effects of
#'     the exposures on the outcome.}
#'   \item{theta.cov}{A (p x p) covariance matrix of \code{theta} obtained
#'     from the bootstrap samples (for multivariable input).}
#'   \item{theta.se}{A vector (length p) of the standard errors of
#'     \code{theta}.}
#'   \item{gamma}{A vector (length n) of the estimated pleiotropy or
#'     background bias terms for each SNP.}
#'   \item{theta.bootstrap}{A (sampling.time x p) matrix of bootstrap
#'     estimates of \code{theta}.}
#' }
#'
#' For univariable input (when \code{bX} is a vector), \code{Grapple_BBC}
#' calls \code{MRRAPS_BBC} internally and returns its output (scalar effect,
#' variance, standard error, pleiotropy term, and bootstrap replicates).
#'
#' @importFrom MASS rlm
#' @import CppMatrix
#' @export

GRAPPLE_BBC=function(by,bX,byse,bXse,cov_RB,max.iter=100,max.eps=1e-6,lambda=3,a=3,tau_upper=10,sampling.time=300,n_threads=4,max.prop.pleio=0.5,sampling.strategy="subsampling"){
if(is.vector(bX)==T){
A=MRRAPS_BBC(by=by,bx=bX,byse=byse,bxse=bXse,cov_RB=cov_RB,max.iter=max.iter,max.eps=max.eps,lambda=lambda,a=a,tau_upper=tau_upper,sampling.time=sampling.time,n_threads=n_threads,max.prop.pleio=max.prop.pleio,sampling.strategy=sampling.strategy)
}else{
######### Basic Processing  ##############
fit.ini=MRBEE_BBC(by=by,bX=bX,byse=byse,bXse=bXse,cov_RB=cov_RB,max.iter=max.iter,max.eps=max.eps,sampling.time=10)
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
eta=MCP_simulation(iter=5e5,lambda=lambda,gamma=a)
n=m=length(by)
p=ncol(bX)
RxyList=array(0,c(p+1,p+1,m))
for(i in 1:m){
A=cov_RB[[i]]*byseinv[i]^2
RxyList[,,i]=A
}
########## Initial Estimation ############
theta.ini=fit.ini$theta
theta=theta.ini
theta1=10000
e=c(by-matrixVectorMultiply(bX,theta))
gamma.ini=soft(e,3)
gamma=gamma.ini
tau=0.1
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-matrixVectorMultiply(bX,theta))
gra_stat=grapple_stat(RxyList=RxyList,theta=theta,e=e-gamma,n_threads=n_threads)
var_vec=as.vector(gra_stat$var_vec)+tau
var_vecinv=1/var_vec
var_cor=gra_stat$var_cor
bias_correction=gra_stat$bias_correction
gamma=mcp(e,lam=lambda*pmax(1,sqrt(var_vec)),a=a)
gamma=pleio_adj(gamma,max.prop.pleio)
g_theta=-matrixVectorMultiply(t(bX),(e-gamma)*var_vec)-matrixVectorMultiply(t(var_cor),(e-gamma)^2*var_vecinv^2)
Hinv=matrixMultiply(t(bX),bX*(var_vecinv))-bias_correction
Hinv=matrixInverse(Hinv)
theta=theta-as.vector(Hinv%*%g_theta)*min(max(1,1.1-iter/10),0.2)
tau_eq <- function(tau) {
z <- (e - gamma) / sqrt(gra_stat$var_vec + tau)
sum(rho_mcp(z, lambda = lambda, gamma = a)) - m * eta
}
tau <- solve_tau(tau_eq, tau_upper)
iter=iter+1
if(iter>10) error=sqrt(sum((theta-theta1)^2))/sqrt(length(theta))
}

################# Inference ##############
ThetaMatrix=matrix(0,sampling.time,p)
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:m,m,replace=T)
mi=m
}else{
ind=sample(1:m,m*0.5,replace=F)
mi=0.5*m
}
thetai=theta*runif(p,0.95,1.05)
theta1i=10000
bXi=bX[ind,]
byi=by[ind]
gammai=gamma[ind]
RxyListi <- RxyList[,,ind]
taui=tau
errori=sqrt(sum((thetai-theta1i)^2))
iteri=0
while(errori>max.eps&iteri<max.iter){
theta1i=thetai
ei=c(byi-matrixVectorMultiply(bXi,thetai))
gra_stati=grapple_stat(RxyList=RxyListi,theta=thetai,e=ei-gammai,n_threads = n_threads)
var_veci=gra_stati$var_vec+taui
var_cori=gra_stati$var_cor
var_vecinvi=1/var_veci
bias_correctioni=gra_stati$bias_correction
gammai=mcp(ei,lam=pmax(1,lambda*sqrt(var_veci)),a=a)
gammai=pleio_adj(gammai,max.prop.pleio)
g_thetai=-matrixVectorMultiply(t(bXi),(ei-gammai)*var_vecinvi)-matrixVectorMultiply(t(var_cori),(ei-gammai)^2*var_vecinvi^2)
Hinvi=matrixMultiply(t(bXi),bXi*var_vecinvi)-bias_correctioni
Hinvi=matrixInverse(Hinvi)
thetai=thetai-as.vector(Hinvi%*%g_thetai)*min(max(1,1.1-iteri/10),0.2)
tau_eqi <- function(taui) {
z <- (ei - gammai) / sqrt(gra_stati$var_vec + taui)
sum(rho_mcp(z, lambda = lambda, gamma = a)) - mi * eta
}
taui <- solve_tau(tau_eqi, tau_upper)
iteri=iteri+1
if(iteri>7) errori=sqrt(sum((thetai-theta1i)^2))/sqrt(length(theta))
}
ThetaMatrix[i,]=thetai
}
theta.cov=cov(ThetaMatrix)
theta.se=sqrt(diag(theta.cov))
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
A=list()
A$theta=c(theta)
A$theta.cov=theta.cov
A$theta.se=theta.se
A$gamma=gamma/byseinv
A$theta.bootstrap=ThetaMatrix
A$tau=tau
}
return(A)
}
