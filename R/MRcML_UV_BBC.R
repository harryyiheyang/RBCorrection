#' Univariable MRcML with background bias correction (BBC).
#'
#' This function estimates the univariable causal effect using an MRcML-type
#' likelihood method with background bias correction (BBC). It uses
#' an MCP-based loss to downweight outlying instruments and incorporates
#' LDSC and Rao–Blackwell corrections to improve stability under weak and
#' biased instruments.
#'
#' @param by A numeric vector (m x 1) of GWAS effect sizes for the outcome
#'   after Rao–Blackwell correction.
#' @param bx A numeric vector (m x 1) of GWAS effect sizes for the exposure
#'   after Rao–Blackwell correction.
#' @param byse A numeric vector (m x 1) of standard errors for \code{by}
#'   after Rao–Blackwell correction.
#' @param bxse A numeric vector (m x 1) of standard errors for \code{bx}
#'   after Rao–Blackwell correction.
#' @param cov_RB A list of length m, where each element is a 2 x 2 matrix of
#'   Rao–Blackwell correction terms for the joint (exposure, outcome)
#'   GWAS estimates at a given SNP.
#' @param max.iter Maximum number of iterations for updating the causal
#'   effect. Default is 30.
#' @param max.eps Tolerance for the stopping criterion based on the change
#'   in the causal effect estimate. Default is 1e-5.
#' @param lambda Tuning parameter for the MCP loss used to downweight
#'   large residuals. Default is 3.
#' @param a Concavity parameter for the MCP loss (often denoted gamma in
#'   the MCP literature). Must be greater than 1. Default is 5.
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   variance of the causal effect. Default is 300.
#' @param n_threads Number of threads used in the internal C++ routines.
#'   Defaults to 4.
#' @param max.prop.pleio A numeric value in (0, 1) giving the maximum allowed
#'   proportion of non-zero pleiotropic terms. Default is 0.5.
#' @param sampling.strategy Character string indicating the resampling
#'   strategy, either \code{"bootstrap"} (sampling with replacement) or
#'   \code{"subsampling"} (using 50\% of the instruments without replacement).
#'
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

MRcML_UV_BBC=function(by,bx,byse,bxse,cov_RB,max.iter=50,max.eps=1e-6,lambda=3,a=3,sampling.time=300,n_threads=4,max.prop.pleio=0.5,sampling.strategy="subsampling"){
######### Basic Processing  ##############
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=m=length(by)
p=1
ThetaList=array(0,c(p+1,p+1,m))
for(i in 1:m){
A=cov_RB[[i]]*byseinv[i]^2
ThetaList[,,i]=solve(A)
}
########## Initial Estimation ############
fit.ini=rlm(by~bx-1)
theta.ini=fit.ini$coefficients
theta=theta.ini
theta1=10000
e=c(by-bx*theta)
gamma.ini=soft(e,3)
gamma=gamma.ini
bxest=bx
########## Iteration ###################
error=abs(theta-theta1)
iter=0
Thetayy=ThetaList[p+1,p+1,]
ThetaX2y=ThetaList[1:p,p+1,]
while(error>max.eps&iter<max.iter){
theta1=theta
e=by-gamma
bxest=MRcML_UV_bxest(ThetaList,bx,e,theta,n_threads)
H=sum(bxest^2*Thetayy)
g=sum(bxest*(by-gamma)*Thetayy)+sum(bxest*(bx-bxest)*ThetaX2y)
theta=g/H
res=(bx-bxest)*ThetaX2y/Thetayy+by-bxest*theta
gamma=mcp(res,lam=lambda/Thetayy,a=a)
gamma=pleio_adj(gamma,max.prop.pleio)
iter=iter+1
if(iter>10) error=abs(theta-theta1)
}
################# Inference ##############
ThetaVec=c(1:sampling.time)*0
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:m,m,replace=T)
}else{
ind=sample(1:m,m*0.5,replace=F)
}
thetai=theta*runif(p,0.95,1.05)
theta1i=10000
bxi=bx[ind]
bxesti=bx[ind]
byi=by[ind]
gammai=gamma[ind]
ThetaListi=ThetaList[,,ind]
errori=abs(thetai-theta1i)
iteri=0
ThetaX2yi=ThetaListi[1:p,p+1,]
Thetayyi=ThetaListi[p+1,p+1,]
while(errori>max.eps&iteri<max.iter){
theta1i=thetai
ei=byi-gammai
bxesti=MRcML_UV_bxest(ThetaListi,bxi,ei,thetai,n_threads)
Hi=sum(bxesti^2*Thetayyi)
gi=sum(bxesti*(byi-gammai)*Thetayyi)+sum(bxesti*(bxi-bxesti)*ThetaX2yi)
thetai=gi/Hi
resi=(bxi-bxesti)*ThetaX2yi/Thetayyi+byi-bxesti*thetai
gammai=mcp(resi,lam=lambda/Thetayyi,a=a)
gammai=pleio_adj(gammai,max.prop.pleio)
iteri=iteri+1
if(iteri>7) errori=abs(thetai-theta1i)
}
ThetaVec[i]=thetai
}
theta.var=var(ThetaVec)
theta.se=sqrt(theta.var)
names(gamma)=names(bx)
A=list()
A$theta=c(theta)
A$theta.var=theta.var
A$theta.se=theta.se
A$gamma=gamma/byseinv
A$theta.bootstrap=ThetaVec
return(A)
}
