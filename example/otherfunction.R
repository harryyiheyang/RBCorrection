
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

imrpdetect=function(x,theta,RxyList,indvalid,var.est="robust",FDR=T,adjust.method="Sidak"){
  p=length(theta)
  if(var.est=="robust"){
    varx=stats::mad(x[indvalid])^2
  }
  if(var.est=="variance"){varx=stats::var(x[indvalid])}
  if(var.est=="ordinal"){
    varx=x*0
    for(i in 1:length(x)){
      varx[i]=c(RxyList[i,p+1,p+1]+t(theta)%*%RxyList[i,1:p,1:p]%*%theta-2*sum(theta*RxyList[i,p+1,1:p]))
    }
  }
  pv=stats::pchisq(x^2/varx,1,lower.tail=F)
  if(FDR==T){
    pv=FDRestimation::p.fdr(pvalues=pv,adjust.method=adjust.method)$fdrs
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

biasterm=function(RxyList,indvalid){
  X=RxyList[1,,]*0
  for(i in indvalid){
    X=X+RxyList[i,,]
  }
  return(X)
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

MRBEE_UV_Winner=function (by,bx,byse,bxse,Rxy,S_RB,max.iter=30,max.eps=1e-04,
                          pv.thres=0.05,phi=1,var.est="variance",FDR=T,
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
S_RB_Sum <- Reduce("+", Map(function(S, w) S * w, S_RB[indvalid], byseinv[indvalid]^2))
h=sum(bx[indvalid]^2)-sum(bxse[indvalid]^2*Rxy[1,1])+S_RB_Sum[1,1]
h=h+exp(phi-mu_min/sqrt(m))/h
g=sum(bx[indvalid]*by[indvalid])-Rxy[1,2]*sum(bxse[indvalid] *byse[indvalid])+S_RB_Sum[1,2]
theta=g/h
iter=iter+1
if (iter>5)
error=sqrt(sum((theta-theta1)^2))
}
S_RBArray=abind::abind(S_RB,  along = 3)
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


MRcML_UV_Winner=function(by,bx,byse,bxse,Rxy,S_RB,max.iter=50,max.eps=1e-04,
                          Kvec=c(0:round(length(bx)/4)),n=length(bx),
                         sampling.iter=15,sampling.time=1000,theta.ini=NULL){
m=length(by)
###############################################################################
get_thetaxy=function(bxse,byse,Rxy,S_RB=NULL){
m=length(bxse)
ThetaList=array(0,c(m,2,2))
Thetaxy=solve(Rxy)
if(is.null(S_RB)){
for(i in 1:m){
zi=c(bxse[i],byse[i])
ziinv=zi
ziinv[ziinv>0]=1/ziinv[ziinv>0]
ThetaList[i,,]=t(t(Thetaxy)*ziinv)*ziinv
}
}else{
for(i in 1:m){
zi=c(bxse[i],byse[i])
Rxyi=t(t(Rxy)*zi)*zi
ThetaList[i,,]=matrixGeneralizedInverse(Rxyi-S_RB[[i]])
}
}
return(ThetaList)
}
###############################################################################
ThetaList=get_thetaxy(bxse,byse,Rxy,S_RB)
if(is.null(theta.ini)){
fit.ini=MASS::rlm(by~bx-1,weights=1/byse^2)
e.ini=abs(fit.ini$residuals/byse)
theta.ini=coef(fit.ini)
}else{
e.ini=by-bx*theta.ini
e.ini=abs(e.ini/byse)
}
bZ=cbind(bx,by)
Bic=Btheta=Kvec
Be=matrix(0,length(Kvec),m)
for(k in 1:length(Kvec)){
e=e.ini
error=1
theta=theta.ini
for(iter in 1:max.iter){
theta_prev=theta
vartheta=c(1,theta)
indvalid=get_top_indices(e,Kvec[k])
################################################################################
fit.iter=update_theta_e(by,bx,byse,ThetaList,vartheta,indvalid)
theta=fit.iter$theta
e=fit.iter$e
if(iter>3){
error=abs(theta-theta_prev)
}
if(error<max.eps) break
}
Bic[k]=log(fit.iter$loss)+k/n*log(n)
Btheta[k]=theta
}
kstar=which.min(Bic)
theta=Btheta[kstar]
e=Be[kstar,]
indvalid=get_top_indices(e,Kvec[kstar])
############################### Inference #######################################
thetavec=c(1:sampling.time)
for(j in 1:sampling.time){
ind=sample(m,m,replace=T)
byj=by[ind]
bxj=bx[ind]
bysej=byse[ind]
bxsej=bxse[ind]
ej=e[ind]
thetaj=theta*runif(1,0.95,1.05)
bZj=cbind(bxj,byj)
ThetaListj=ThetaList[ind,,]
errorj=1
for(iterj in 1:sampling.iter){
theta_prevj=thetaj
varthetaj=c(1,thetaj)
indvalidj=get_top_indices(ej,Kvec[kstar])
fit.iterj=update_theta_e(byj,bxj,bysej,ThetaListj,varthetaj,indvalidj)
thetaj=fit.iterj$theta
ej=fit.iterj$e
if(iterj>3){
errorj=abs(thetaj-theta_prevj)
}
if(errorj<max.eps) break
}
thetavec[j]=thetaj
}
theta.se=sd(thetavec)
A=list(theta=theta,theta.se=theta.se,theta.vec=thetavec,k.opt=kstar,BICvec=Bic,IV.valid=indvalid)
return(A)
}
