#' MRcML with background bias correction (BBC).
#'
#' This function estimates multivariable causal effects using an MRcML-type
#' likelihood method combined with background bias correction (BBC).
#' It uses an MCP-based loss to downweight outlying instruments and
#' accounts for sample overlap through LDSC and Rao–Blackwell corrections
#' applied to both exposure and outcome GWAS estimates.
#'
#' @param by A vector (n x 1) of the GWAS effect sizes of the outcome
#'   after Rao–Blackwell correction.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures
#'   after Rao–Blackwell correction. If univariable, \code{bX} can be a
#'   vector of length n, in which case \code{MRcML_BBC} calls
#'   \code{MRcML_UV_BBC} internally.
#' @param byse A vector (n x 1) of the standard errors of \code{by}
#'   after Rao–Blackwell correction.
#' @param bXse A matrix (n x p) of the standard errors of \code{bX}
#'   after Rao–Blackwell correction. If univariable, \code{bXse} can be
#'   a vector of length n.
#' @param gcov A ((p+1) x (p+1)) matrix of the per-SNP genetic covariance
#'   of the p exposures and the outcome. The last row/column corresponds
#'   to the outcome.
#' @param ldsc A vector (n x 1) of LD scores of the instruments.
#' @param cov_RB A list of length n, where each element is a ((p+1) x (p+1))
#'   matrix of Rao–Blackwell correction terms for the joint
#'   (exposures, outcome) GWAS estimates at a given SNP.
#' @param max.iter Maximum number of iterations for causal effect estimation.
#'   Defaults to 30.
#' @param max.eps Tolerance for the stopping criterion based on the change
#'   in the causal effect estimates. Defaults to 1e-5.
#' @param lambda Tuning parameter for the MCP loss used to downweight
#'   outlying instruments. Defaults to 3.
#' @param a Concavity parameter of the MCP loss (often denoted gamma in the
#'   MCP literature). Must be greater than 1. Defaults to 5.
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   covariance of the causal effect estimates. Defaults to 300.
#' @param n_threads Number of threads used in the internal C++ routines.
#'   Defaults to 4.
#' @param max.prop.pleio A numeric value in (0, 1) giving the maximum allowed
#'   proportion of non-zero pleiotropic terms. Instruments beyond this
#'   proportion are truncated by \code{pleio_adj}. Default is 0.5.
#' @param sampling.strategy Character string indicating the resampling
#'   strategy, either \code{"bootstrap"} (sampling with replacement) or
#'   \code{"subsampling"} (using 50\% of the instruments without replacement).
#'
#' @return A list containing:
#' \describe{
#'   \item{theta}{A vector (length p) of the estimated causal effects of
#'     the exposures on the outcome.}
#'   \item{theta.cov}{A (p x p) covariance matrix of \code{theta} obtained
#'     from the bootstrap samples.}
#'   \item{theta.se}{A vector (length p) of the standard errors of
#'     \code{theta}.}
#'   \item{gamma}{A vector (length n) of the estimated pleiotropic or
#'     background bias terms for each SNP.}
#'   \item{theta.bootstrap}{A (sampling.time x p) matrix of bootstrap
#'     estimates of \code{theta}.}
#' }
#'
#' For univariable input (when \code{bX} is a vector), \code{MRcML_BBC}
#' calls \code{MRcML_UV_BBC} internally and returns its output.
#'
#' @importFrom MASS rlm
#' @import CppMatrix
#' @export

MRcML_BBC=function(by,bX,byse,bXse,gcov,ldsc,cov_RB,max.iter=30,max.eps=1e-5,lambda=3,a=3,sampling.time=300,n_threads=4,max.prop.pleio=0.5,sampling.strategy="subsampling"){
if(is.vector(bX)==T){
A=MRcML_UV_BBC(by=by,bx=bX,byse=byse,bxse=bXse,ldsc=ldsc,cov_RB=cov_RB,max.iter=max.iter,max.eps=max.eps,lambda=lambda,a=a,sampling.time=sampling.time,n_threads=n_threads,max.prop.pleio=max.prop.pleio,sampling.strategy=sampling.strategy)
}else{
######### Basic Processing  ##############
fit.ini=MRBEE_BBC(by=by,bX=bX,byse=byse,bXse=bXse,cov_RB=cov_RB,gcov=gcov,ldsc=ldsc,max.iter=max.iter,max.eps=max.eps)
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=m=length(by)
p=ncol(bX)
ThetaList=array(0,c(p+1,p+1,m))
for(i in 1:m){
A=(cov_RB[[i]]+ldsc[i]*gcov)*byseinv[i]^2
ThetaList[,,i]=solve(A)
}
########## Initial Estimation ############
theta.ini=fit.ini$theta
theta=theta.ini
theta1=10000
e=c(by-matrixVectorMultiply(bX,theta))
gamma.ini=soft(e,3)
gamma=gamma.ini
bXest=bX
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
Thetayy=ThetaList[p+1,p+1,]
ThetaX2y=ThetaList[1:p,p+1,]
while(error>max.eps&iter<max.iter){
theta1=theta
e=by-gamma
bXest=MRcML_bXest(ThetaList,bX,e,theta,n_threads)
H=matrixMultiply(t(bXest),bXest*Thetayy)
g=matrixVectorMultiply(t(bXest),(by-gamma)*Thetayy+rowSums((bX-bXest)*t(ThetaX2y)))
theta=as.vector(solve(H)%*%g)
res=rowSums((bX-bXest)*t(ThetaX2y))/Thetayy+by-matrixVectorMultiply(bXest,theta)
gamma=mcp(res,lam=lambda/Thetayy*byse,a=a)
gamma=pleio_adj(gamma,max.prop.pleio)
iter=iter+1
if(iter>10) error=sqrt(sum((theta-theta1)^2))
}
################# Inference ##############
ThetaMatrix=matrix(0,sampling.time,p)
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:m,m,replace=T)
}else{
ind=sample(1:m,m*0.5,replace=F)
}
thetai=theta*runif(p,0.95,1.05)
theta1i=10000
bXi=bX[ind,]
bXesti=bX[ind,]
byi=by[ind]
gammai=gamma[ind]
ThetaListi=ThetaList[,,ind]
errori=sqrt(sum((thetai-theta1i)^2))
iteri=0
ThetaX2yi=ThetaListi[1:p,p+1,]
Thetayyi=ThetaListi[p+1,p+1,]
while(errori>max.eps&iteri<max.iter){
theta1i=thetai
ei=byi-gammai
bXesti=MRcML_bXest(ThetaListi,bXi,ei,thetai,n_threads)
Hi=matrixMultiply(t(bXesti),bXesti*Thetayyi)
gi=matrixVectorMultiply(t(bXesti),(byi-gammai)*Thetayyi+rowSums((bXi-bXesti)*t(ThetaX2yi)))
thetai=as.vector(solve(Hi)%*%gi)
resi=rowSums((bXi-bXesti)*t(ThetaX2yi))/Thetayyi+byi-matrixVectorMultiply(bXesti,thetai)
gammai=mcp(resi,lam=lambda/Thetayyi*byse[ind],a=a)
gammai=pleio_adj(gammai,max.prop.pleio)
iteri=iteri+1
if(iteri>7) errori=sqrt(sum((thetai-theta1i)^2))
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
}
return(A)
}

