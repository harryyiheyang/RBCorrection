#' MRBEE with background bias correction (BBC).
#'
#' This function estimates the causal effect using a bias-correction estimating equation, combining with background bias correction (BBC).
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome yielded by Rao-Blackwell Correction.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures yielded by Rao-Blackwell Correction.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome yielded by Rao-Blackwell Correction.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposuresyielded by Rao-Blackwell Correction.
#' @param cov_RB A list of m matrices of Rao-Blackwell correction terms, each of dimension (p+1) × (p+1).
#' @param gcov A matrix (p+1 x p+1) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 30.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param pv.thres P-value threshold in pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the variance of residual in pleiotropy test. Can be "robust", "variance", or "ordinal". Defaults is "variance" that estimates the variance of residual using median absolute deviation (MAD).
#' @param FDR Logical. Whether to apply the FDR to convert the p-value to q-value. Defaults to TRUE.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   covariance of the causal effect estimates. Defaults to 300.
#' @param sampling.strategy Character string indicating the resampling
#'   strategy, either \code{"bootstrap"} (sampling with replacement) or
#'   \code{"subsampling"} (using 50\% of the instruments without replacement).
#' @param n_threads The number of threads. Defaults to 4.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy
#' @importFrom MASS rlm
#' @import CppMatrix
#' @importFrom abind abind
#' @export
#'
MRBEE_BBC=function(by,bX,byse,bXse,gcov,ldsc,cov_RB,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak",sampling.time=300,sampling.strategy="subsampling",n_threads=4){
if(is.vector(bX)==T){
A=MRBEE_UV_BBC(by=by,bx=bX,byse=byse,bxse=bXse,ldsc=ldsc,max.iter=max.iter,max.eps=max.eps,pv.thres=pv.thres,var.est=var.est,FDR=FDR,adjust.method=adjust.method,cov_RB=cov_RB,sampling.time=sampling.time,sampling.strategy=sampling.strategy,n_threads=n_threads)
A$gamma=A$delta
A$theta.se=sqrt(A$vartheta)
}else{
######### Basic Processing  ##############
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=m=length(by)
p=ncol(bX)
RxyList=array(0,c(p+1,p+1,m))
for(i in 1:m){
A=(cov_RB[[i]]+ldsc[i]*gcov)*byseinv[i]^2
RxyList[,,i]=A
}
Rxyall=biasterm(RxyList=RxyList,c(1:n),n_threads=n_threads)
########## Initial Estimation ############
fit=MASS::rlm(by~bX-1)
theta.ini=fit$coefficient
theta=theta.ini
theta1=10000
e=c(by-matrixVectorMultiply(bX,theta))
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5) ## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-matrixVectorMultiply(bX,theta))
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>=pv.thres)
if (length(indvalid) <= length(pv) * 0.5) {
indvalid.cut = which(pv >= stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid),n_threads)
}
Hinv=matrixMultiply(t(bX[indvalid,]),bX[indvalid,])-Rxysum[1:p,1:p]
Hinv=matrixInverse(Hinv)
g=matrixVectorMultiply(t(bX[indvalid,]),by[indvalid])-Rxysum[1+p,1:p]
theta=matrixVectorMultiply(Hinv,g)
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
e[indvalid]=0
################# Inference ##############
ThetaMatrix=matrix(0,sampling.time,p)
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:n,n,replace=T)
mi=n
}else{
ind=sample(1:n,n*0.5,replace=F)
mi=n/2
}
thetai=theta*runif(p,0.95,1.05)
theta1i=10000
bXi=bX[ind,]
byi=by[ind]
indvalidi=which(e[ind]==0)
RxyListi=RxyList[,,ind]
errori=sqrt(sum((thetai-theta1i)^2))
iteri=0
while(errori>max.eps&iteri<max.iter){
theta1i=thetai
ei=c(byi-matrixVectorMultiply(bXi,thetai))
pvi=imrpdetect(x=ei,theta=thetai,RxyList=RxyListi,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalidi)
indvalidi=which(pvi>=pv.thres)
if (length(indvalidi) <= length(pvi) * 0.5) {
indvalid.cut = which(pvi >= stats::quantile(pvi, 0.5))
indvalidi = union(indvalidi, indvalid.cut)
}
Rxysumi=biasterm(RxyList=RxyListi,indvalidi,n_threads)
Hinvi=matrixMultiply(t(bXi[indvalidi,]),bXi[indvalidi,])-Rxysumi[1:p,1:p]
Hinvi=matrixInverse(Hinvi)
gi=matrixVectorMultiply(t(bXi[indvalidi,]),byi[indvalidi])-Rxysumi[1+p,1:p]
thetai=matrixVectorMultiply(Hinvi,gi)
iteri=iteri+1
if(iteri>5) errori=sqrt(sum((thetai-theta1i)^2))
}
ThetaMatrix[i,]=thetai
}
theta.cov=cov(ThetaMatrix)
theta.se=sqrt(diag(theta.cov))
names(theta)=colnames(bX)
r=(by-matrixVectorMultiply(bX,theta))*byse1
r[indvalid]=0
names(theta)=colnames(bX)
names(r)=rownames(bX)
A=list()
A$theta=c(theta)
A$theta.cov=theta.cov
A$theta.se=theta.se
A$gamma=r
}
return(A)
}
