#' Univariable MRBEE with background bias correction (BBC)
#'
#' This function estimates the univariable causal effect using bias-corrected estimating equations with background bias correction (BBC) to improve numerical stability under weak instruments.
#'
#' @param by A numeric vector (m × 1) of GWAS effect sizes for the outcome yielded by Rao-Blackwell Correction.
#' @param bx A numeric vector (m × 1) of GWAS effect sizes for the exposure yielded by Rao-Blackwell Correction.
#' @param byse A numeric vector (m × 1) of standard errors for `by` yielded by Rao-Blackwell Correction.
#' @param bxse A numeric vector (m × 1) of standard errors for `bx` yielded by Rao-Blackwell Correction.
#' @param cov_RB A list of m matrices of Rao-Blackwell correction terms, each of dimension 2 × 2.
#' @param max.iter Maximum number of iterations for updating the causal effect. Default is 30.
#' @param max.eps Tolerance for stopping criteria. Default is 1e-4.
#' @param pv.thres P-value threshold used in pleiotropy detection. Default is 0.05.
#' @param var.est Method for estimating the variance of residuals. One of "robust", "variance", or "ordinal". Default is "variance".
#' @param FDR Logical. Whether to convert p-values to q-values using FDR adjustment. Default is TRUE.
#' @param adjust.method Method used to adjust p-values. Default is "Sidak".
#' @param sampling.time Number of bootstrap iterations used to estimate the
#'   variance of the causal effect. Default is 300.
#' @param n_threads The number of threads. Defaults to 4.
#' @param sampling.strategy "bootstrap" or "subsampling" (0.5 sample without replacement).

#' @return A list with the following components:
#' \item{theta}{Estimated causal effect.}
#' \item{vartheta}{Estimated variance of theta.}
#' \item{delta}{Estimated pleiotropic residuals (set to zero for valid IVs).}
#' @importFrom MASS rlm
#' @importFrom abind abind
#' @export

MRBEE_UV_BBC=function(by,bx,byse,bxse,cov_RB,max.iter=50,max.eps=1e-06,sampling.time=300,sampling.strategy="subsampling",n_threads,
                    pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak"){
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
RxyList=array(0,c(2,2,m))
for(i in 1:m){
RxyList[,,i]=cov_RB[[i]]*byseinv[i]^2
}
error=abs(theta-theta1)
iter=0
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
var_e=mean(e[indvalid]^2)
h=sum(bx[indvalid]^2)-sum(RxyList[1,1,indvalid]*e[indvalid]^2/var_e)
g=sum(bx[indvalid]*by[indvalid])-sum(RxyList[1,2,indvalid]*e[indvalid]^2/var_e)
theta=g/h
iter=iter+1
if (iter>8) error=abs(theta-theta1)
}
e[indvalid]=0
################# Inference ##############
ThetaVec=c(1:sampling.time)*0
for(i in 1:sampling.time){
if(sampling.strategy=="bootstrap"){
ind=sample(1:m,m,replace=T)
}else{
ind=sample(1:m,m*0.5,replace=F)
}
thetai=theta*runif(1,0.95,1.05)
theta1i=10000
bxi=bx[ind]
byi=by[ind]
indvalidi=which(e[ind]==0)
RxyListi=RxyList[,,ind]
errori=abs(thetai-theta1i)
iteri=0
while (errori>max.eps & iteri<max.iter) {
theta1i=thetai
ei=c(byi-bxi*thetai)
pvi=imrpdetect(x=ei,theta=thetai,RxyList=RxyListi,
              var.est=var.est,FDR=FDR,adjust.method=adjust.method,
              indvalid=indvalidi)
pvi[is.na(pvi)]=1
indvalidi=which(pvi>=pv.thres)
if (length(indvalidi)<=length(pvi)*0.5) {
  indvalid.cut=which(pvi>=stats::quantile(pvi,0.5))
  indvalidi=union(indvalidi,indvalid.cut)
}
var_ei=mean(ei[indvalidi]^2)
hi=sum(bxi[indvalidi]^2)-sum(RxyListi[1,1,indvalidi]*ei[indvalidi]^2/var_ei)
gi=sum(bxi[indvalidi]*byi[indvalidi])-sum(RxyListi[1,2,indvalidi]*e[indvalidi]^2/var_ei)
thetai=gi/hi
iteri=iteri+1
if (iteri>7) errori=abs(thetai-theta1i)
}
ThetaVec[i]=thetai
}
theta.var=var(ThetaVec)
theta.se=sqrt(theta.var)

A=list()
A$theta=theta
A$theta.var=theta.var
A$theta.se=theta.se
r=c(by-bx*theta)*byse1
r[indvalid]=0
names(r)=rownames(bx)
A$delta=r
return(A)
}
