#' Univariable MRBEE with background bias correction (BBC)
#'
#' This function estimates the univariable causal effect using bias-corrected estimating equations with background bias correction (BBC) to improve numerical stability under weak instruments.
#'
#' @param by A numeric vector (m × 1) of GWAS effect sizes for the outcome yielded by Rao-Blackwell Correction.
#' @param bx A numeric vector (m × 1) of GWAS effect sizes for the exposure yielded by Rao-Blackwell Correction.
#' @param byse A numeric vector (m × 1) of standard errors for `by` yielded by Rao-Blackwell Correction.
#' @param bxse A numeric vector (m × 1) of standard errors for `bx` yielded by Rao-Blackwell Correction.
#' @param cov_RB A list of m matrices of Rao-Blackwell correction terms, each of dimension 2 × 2.
#' @param gcov A matrix (2 x 2) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.
#' @param max.iter Maximum number of iterations for updating the causal effect. Default is 30.
#' @param max.eps Tolerance for stopping criteria. Default is 1e-4.
#' @param pv.thres P-value threshold used in pleiotropy detection. Default is 0.05.
#' @param phi Spectrum regularization parameter. Controls the shrinkage strength. Default is 2.
#' @param var.est Method for estimating the variance of residuals. One of "robust", "variance", or "ordinal". Default is "variance".
#' @param FDR Logical. Whether to convert p-values to q-values using FDR adjustment. Default is TRUE.
#' @param adjust.method Method used to adjust p-values. Default is "Sidak".
#'
#' @return A list with the following components:
#' \item{theta}{Estimated causal effect.}
#' \item{vartheta}{Estimated variance of theta.}
#' \item{delta}{Estimated pleiotropic residuals (set to zero for valid IVs).}
#' @importFrom MASS rlm
#' @importFrom abind abind
#' @export

MRBEE_UV_BBC=function(by,bx,byse,bxse,cov_RB,gcov=diag(2)*0,ldsc=rep(0,length(by)),max.iter=30,max.eps=1e-04,
                      pv.thres=0.05,phi=2,var.est="variance",FDR=T,adjust.method="Sidak"){
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=m=length(by)
fit=MASS::rlm(by~bx-1)
theta=fit$coefficient
theta1=10000
e=c(by-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
RxyList=array(0,c(n,2,2))
for(i in 1:m){
RxyList[i,,]=(cov_RB[[i]]+ldsc[i]*gcov)*byseinv[i]^2
}
error=abs(theta-theta1)
iter=0
mu_min=max(0,sum(bx^2/RxyList[,1,1])-m)
while (error>max.eps & iter<max.iter) {
theta1=theta
e=c(by-bx*theta)
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,
              var.est=var.est,FDR=FDR,adjust.method=adjust.method,
              indvalid=indvalid)
pv[is.na(pv)]=1
indvalid=which(pv>=pv.thres)
if (length(indvalid)<=length(pv)*0.5) {
indvalid.cut=which(pv>=stats::quantile(pv,0.5))
indvalid=union(indvalid,indvalid.cut)
}
h=sum(bx[indvalid]^2)-sum(RxyList[indvalid,1,1])
h=h+exp(phi-mu_min/sqrt(m))/h
g=sum(bx[indvalid]*by[indvalid])-sum(RxyList[indvalid,1,2])
theta=g/h
iter=iter+1
if (iter>5) error=abs(theta-theta1)
}
adjf=n/(length(indvalid)-1)
Hat=outer(bx[indvalid],bx[indvalid])/h
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
e[indvalid]=e[indvalid]/Hat
E = -bx[indvalid]*e[indvalid] + RxyList[indvalid,1,2] - RxyList[indvalid,1,1]*theta
vartheta=sum(E^2)/h^2*adjf

A=list()
A$theta=theta
A$vartheta=vartheta
A$theta.se=sqrt(vartheta)
r=c(by-bx*theta)*byse1
r[indvalid]=0
names(r)=rownames(bx)
A$delta=r
return(A)
}
