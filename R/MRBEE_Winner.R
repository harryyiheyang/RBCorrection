#' Mendelian randomization with bias-correction estimating equation: detecting horizontal pleiotropy via hypothesis test.
#'
#' This function estimates the causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The last one should be the outcome.
#' @param S_RB A list of m matrices of Rao-Blackwell correction terms, each of dimension (p+1) × (p+1).
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 30.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param pv.thres P-value threshold in pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the variance of residual in pleiotropy test. Can be "robust", "variance", or "ordinal". Defaults is "variance" that estimates the variance of residual using median absolute deviation (MAD).
#' @param FDR Logical. Whether to apply the FDR to convert the p-value to q-value. Defaults to TRUE.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to 1.5.
#' @param phi The spectrum shrinkage parameter introduced by SRIVW. Default to 2.

#' @return A list containing the estimated causal effect, its covariance, and pleiotropy
#' @importFrom MASS rlm
#' @import CppMatrix
#' @importFrom abind abind
#' @export
#'
MRBEE_Winner=function(by,bX,byse,bXse,Rxy,S_RB,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak",maxdiff=1.5,phi=2){
if(is.vector(bX)==T){
A=MRBEE_UV_Winner(by=by,bx=bX,byse=byse,bxse=bXse,Rxy=Rxy,max.iter=max.iter,max.eps=max.eps,pv.thres=pv.thres,var.est=var.est,FDR=FDR,adjust.method=adjust.method,S_RB=S_SB)
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
n=length(by)
p=ncol(bX)
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
Rxysqrtinv=matrixsqrt(Rxy[1:p,1:p])$wi
G=matrixMultiply(t(bX/bXse),bX/bXse)
G=matrixListProduct(list(Rxysqrtinv,G,Rxysqrtinv))-n*diag(ncol(bX))
mu_min=max(0,min(eigen(G)$values))
S_RBArray=abind::abind(S_RB,  along = 3)
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
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
S_RB_Sum <- Reduce(`+`, lapply(seq_along(indvalid), function(i) {
  S_RBArray[,,indvalid[i]] * byseinv[indvalid[i]]^2
}))
Hinv=matrixMultiply(t(bX[indvalid,]),bX[indvalid,])-Rxysum[1:p,1:p]+S_RB_Sum[1:p,1:p]
Hinv=matrixInverse(Hinv)
g=matrixVectorMultiply(t(bX[indvalid,]),by[indvalid])-Rxysum[1+p,1:p]+S_RB_Sum[1+p,1:p]
theta=matrixVectorMultiply(Hinv,g)
##### MRBEE may generate large theta if Hinv is not well-conditioned. Setting a upper boundary of it.
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
  theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}

iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}

################# Inference ##############
adjf=n/(length(indvalid)-dim(bX)[2])
D=matrixListProduct(list(bX[indvalid,],Hinv,t(bX[indvalid,])))
D=(rep(1,length(indvalid))-diag(D))
D[which(D<0.25)]=0.25
E=-bX[indvalid,]*(e[indvalid]/D)

for(ii in 1:length(indvalid)){
E[ii,]=E[ii,]+RxyList[indvalid[ii],p+1,1:p]-matrixVectorMultiply(RxyList[indvalid[ii],1:p,1:p]-S_RBArray[1:p,1:p, indvalid[ii]],theta)-S_RBArray[1:p,1+p, indvalid[ii]]
}

V=matrixMultiply(t(E),E)
covtheta=matrixListProduct(list(Hinv,V,Hinv))*adjf
r=(by-matrixVectorMultiply(bX,theta))*byse1
r[indvalid]=0
names(theta)=colnames(bX)
names(r)=rownames(bX)
A=list()
A$theta=c(theta)
A$theta.cov=covtheta
A$theta.se=sqrt(diag(covtheta))
A$gamma=r
}
return(A)
}
