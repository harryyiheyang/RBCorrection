#' Univariable MRBEE with Spectrum Regularization
#'
#' This function estimates the univariable causal effect using bias-corrected estimating equations with optional spectrum regularization to improve numerical stability under weak instruments.
#'
#' @param by A numeric vector (m × 1) of GWAS effect sizes for the outcome.
#' @param bx A numeric vector (m × 1) of GWAS effect sizes for the exposure.
#' @param byse A numeric vector (m × 1) of standard errors for `by`.
#' @param bxse A numeric vector (m × 1) of standard errors for `bx`.
#' @param Rxy A 2 × 2 correlation matrix of the exposure and outcome GWAS estimates.
#' @param S_RB A list of m matrices of Rao-Blackwell correction terms, each of dimension 2 × 2.
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

MRBEE_UV_Winner=function (by,bx,byse,bxse,Rxy,S_RB,max.iter=30,max.eps=1e-04,
                    pv.thres=0.05,phi=2,var.est="variance",FDR=T,
                    adjust.method="Sidak"){
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
RxyList=IVweight(byse,bxse,Rxy)
fit=MASS::rlm(by~bx-1)
theta=fit$coefficient
theta1=10000
e=c(by-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
error=abs(theta-theta1)
iter=0
m=length(by)
mu_min=max(0,sum((bx/bxse)^2)/Rxy[1,1]-m)
S_RBArray=abind::abind(S_RB,  along = 3)
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
S_RB_Sum <- Reduce(`+`, lapply(seq_along(indvalid), function(i) {
  S_RBArray[,,indvalid[i]] * byseinv[indvalid[i]]^2
}))
h=sum(bx[indvalid]^2)-sum(bxse[indvalid]^2*Rxy[1,1])+S_RB_Sum[1,1]
h=h+exp(phi-mu_min/sqrt(m))/h
g=sum(bx[indvalid]*by[indvalid])-Rxy[1,2]*sum(bxse[indvalid] *byse[indvalid])+S_RB_Sum[1,2]
theta=g/h
iter=iter+1
if (iter>5)
error=sqrt(sum((theta-theta1)^2))
}
adjf=n/(length(indvalid)-1)
Hat=outer(bx[indvalid],bx[indvalid])/h
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
e[indvalid]=e[indvalid]/Hat
E = -bx[indvalid]*e[indvalid] +
bxse[indvalid]*byse[indvalid] -
bxse[indvalid]^2 * theta +
S_RBArray[1, 1, indvalid] * theta * byseinv[indvalid]^2 -
S_RBArray[1, 2, indvalid] * byseinv[indvalid]^2
vartheta=sum(E^2)/h^2*adjf

A=list()
A$theta=theta
A$vartheta=vartheta
r=c(by-bx*theta)*byse1
r[indvalid]=0
names(r)=rownames(bx)
A$delta=r
return(A)
}
