#' Univariable GRAPPLE with background bias correction (BBC)
#'
#' This function estimates the univariable causal effect using a GRAPPLE-type
#' likelihood method with background bias correction (BBC). It uses
#' an MCP-based loss to downweight outlying instruments and incorporates
#' LDSC and Rao-Blackwell corrections to improve stability under weak and
#' biased instruments.
#'
#' @param by A numeric vector (m x 1) of GWAS effect sizes for the outcome
#'   yielded by Rao-Blackwell correction.
#' @param bx A numeric vector (m x 1) of GWAS effect sizes for the exposure
#'   yielded by Rao-Blackwell correction.
#' @param byse A numeric vector (m x 1) of standard errors for \code{by}
#'   yielded by Rao-Blackwell correction.
#' @param bxse A numeric vector (m x 1) of standard errors for \code{bx}
#'   yielded by Rao-Blackwell correction.
#' @param ldsc A numeric vector (m x 1) of LD scores of the instruments.
#' @param max.iter Maximum number of iterations for updating the causal
#'   effect. Default is 30.
#' @param max.eps Tolerance for the stopping criterion based on the change
#'   in the causal effect estimate. Default is 1e-4.
#' @param lambda Tuning parameter for the MCP loss used to downweight
#'   large residuals. Default is 2.
#' @param a Concavity parameter for the MCP loss (often denoted gamma in
#'   the MCP literature). Must be greater than 1. Default is 3.
#' @param tau_upper Upper bound of the search range for the variance
#'   inflation parameter \code{tau} in the robust variance equation.
#'   Default is 10.
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   variance of the causal effect. Default is 300.
#' @param n_threads The number of threads. Defaults to 4.
#' @param max.prop.pleio A numeric value in (0, 1). The maximum allowed proportion  of non-zero pleiotropic terms. Default is 0.5.
#' @param sampling.strategy "bootstrap" or "subsampling" (0.5 sample without replacement).

#' @return A list with the following components:
#' \item{theta}{Estimated causal effect.}
#' \item{theta.var}{Estimated variance of \code{theta} based on the bootstrap.}
#' \item{theta.se}{Standard error of \code{theta}, computed as
#'   \code{sqrt(theta.var)}.}
#' \item{gamma}{Estimated pleiotropic or background bias term for each SNP.}
#' \item{theta.bootstrap}{A numeric vector of length \code{sampling.time}
#'   containing bootstrap estimates of \code{theta}.}
#'
#' @importFrom MASS rlm
#' @export

MRRAPS_BBC=function(by,bx,byse,bxse,cov_RB,max.iter=100,max.eps=1e-6,lambda=3,a=3,tau_upper=10,sampling.time=300,n_threads=4,max.prop.pleio=0.5,sampling.strategy="subsampling"){
eta = MCP_simulation(iter=5e4, lambda=lambda, gamma=a)
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=m=length(by)
p=1
fit=MASS::rlm(by~bx-1)
theta=fit$coefficient
theta1=10000
e=c(by-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
RxyList=array(0,c(p+1,p+1,m))
for(i in 1:m){
A=cov_RB[[i]]*byseinv[i]^2
RxyList[,,i]=A
}
########## Iteration ###################
gamma.ini=soft(e,1)
gamma.ini[indvalid]=0
tau=0.1
gamma=gamma.ini
error=abs(theta-theta1)
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bx*theta)
gra_stat=grapple_stat(RxyList=RxyList,theta=theta,e=e-gamma,n_threads = n_threads)
var_vec=c(gra_stat$var_vec+tau)
var_cor=c(gra_stat$var_cor)
bias_correction=c(gra_stat$bias_correction)
gamma=mcp(e,lam=lambda*pmax(1,sqrt(var_vec)),a=a)
gamma=pleio_adj(gamma,max.prop.pleio)
g_theta=-sum(bx*(e-gamma)/var_vec)-sum(var_cor*(e-gamma)^2/var_vec^2)
Hinv=sum(bx*bx/var_vec)-bias_correction
theta=theta-g_theta/Hinv*runif(1,0.2,1)
tau_eq <- function(tau) {
z <- (e - gamma) / sqrt(gra_stat$var_vec + tau)
sum(rho_mcp(z, lambda = lambda, gamma = a)) - m * eta
}
tau <- solve_tau(tau_eq, tau_upper)
iter=iter+1
if(iter>10) error=abs(theta-theta1)
}

################# Inference ##############
ThetaVec=c(1:sampling.time)*0
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:m,m,replace=T)
mi=m
}else{
ind=sample(1:m,m*0.5,replace=F)
mi=0.5*m
}
thetai=theta*runif(1,0.95,1.05)
theta1i=10000
bxi=bx[ind]
byi=by[ind]
gammai=gamma[ind]
RxyListi <- RxyList[,,ind]
taui=tau
errori=abs(thetai-theta1i)
iteri=0
while(errori>max.eps&iteri<max.iter){
theta1i=thetai
ei=c(byi-bxi*thetai)
gra_stati=grapple_stat(RxyList=RxyListi,theta=thetai,e=ei-gammai,n_threads = n_threads)
var_veci=c(gra_stati$var_vec+taui)
var_cori=c(gra_stati$var_cor)
bias_correctioni=c(gra_stati$bias_correction)
gammai=mcp(ei,lam=lambda*pmax(1,sqrt(var_veci)),a=a)
gammai=pleio_adj(gammai,max.prop.pleio)
g_thetai=-sum(bxi*(ei-gammai)/var_veci)-sum(var_cori*(ei-gammai)^2/var_veci^2)
Hinvi=sum(bxi*bxi/var_veci)-bias_correctioni
thetai=thetai-g_thetai/Hinvi*runif(1,0.2,1)
tau_eqi <- function(taui) {
z <- (ei - gammai) / sqrt(gra_stati$var_vec + taui)
sum(rho_mcp(z, lambda = lambda, gamma = a)) - mi * eta
}
taui <- solve_tau(tau_eqi, tau_upper)
iteri=iteri+1
if(iteri>10) errori=abs(thetai-theta1i)
}
ThetaVec[i]=thetai
}
theta.var=var(ThetaVec)
theta.se=sqrt(theta.var)
names(gamma)=names(bx)
A=list()
A$theta=theta
A$theta.var=theta.var
A$theta.se=theta.se
A$gamma=gamma/byseinv
A$theta.bootstrap=ThetaVec
A$tau=tau
return(A)
}
