library(MRBEE)
library(MASS)
library(devtools)
library(MRAPSS)
library(dplyr)
library(CppMatrix)
document()
options(bitmapType = "cairo")
m=1000
n=5e4
p=4
Rbb=matrix(0.5,4,4)+diag(4)*0.5
Rxy=matrix(0.5,5,5)+diag(5)*0.5
Rxy=diag(5)
Thetaxx=solve(Rxy[1:p,1:p])
BETA=array(0,c(300,4,5))
theta0=c(0,0.2,0,0)
iv.thes=5e-3
h2=0.05
eta=0.5
i=1

while(i<=300){
indicator=F
tryCatch({
betaX=mvrnorm(n=m,mu=rep(0,p),Rbb)
for(ii in 1:p){
betax=betaX[,ii]
ind=sample(m,0.5*m)
betax[ind]=0
betax=betax/sqrt(sum(betax^2))*sqrt(h2)
betaX[,ii]=betax
}
betay=matrixVectorMultiply(betaX,theta0)+rbinom(m,1,0.1)*rnorm(m,0,sqrt(0.001/0.1))
E=mvrnorm(n=m,mu=rep(0,p+1),Sigma=Rxy)
E=E%*%diag(sqrt(rep(1/n,p+1)))
bX=betaX+E[,1:p]
by=betay+E[,p+1]
bxse=byse=rep(sqrt(1/n),m)
bXse=cbind(bxse,bxse,bxse,bxse)
pv=byse
for(ii in 1:m){
z=bX[ii,]/bXse[ii,]
chisq=sum(z*(Thetaxx%*%z))
pv[ii]=pchisq(chisq,p,lower.tail=F)
}
indselect=which(pv<iv.thes)
fit0=MRBEE.IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,pv.thres=0.05)
fit1=MRBEE.IMRP(by=by[indselect],bX=bX[indselect,],byse=byse[indselect],bXse=bXse[indselect,],Rxy=Rxy,pv.thres=0.05)
BETAMatrix=cbind(bX,by)
SEMatrix=cbind(bXse,byse)
rownames(BETAMatrix)=rownames(SEMatrix)=paste0("V",1:m)
colnames(BETAMatrix)=colnames(SEMatrix)=c(paste0("Exposure",1:p),"Outcome")
###############################################################################
for(ii in 1:m){
z=bX[ii,]/bXse[ii,]+MASS::mvrnorm(n=1,mu=rep(0,p),Sigma=Rxy[1:p,1:p])*eta
chisq=sum(z*(Thetaxx%*%z))
pv[ii]=pchisq(chisq,p,lower.tail=F)
}
indselect=which(pv<iv.thes)
RB=RaoBlackwellCorrect(BETA_Select=BETAMatrix[indselect,],SE_Select=SEMatrix[indselect,],Rxy=Rxy,pv.threshold=iv.thes,eta=eta,B=5000)
fit2=MRBEE_BBC(bX=RB$bX_RB,bXse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,pv.thres=0.05,cov_RB=RB$COV_RB,gcov=diag(p+1)*0,ldsc=rep(0,length(RB$by_RB)))
fit3=GRAPPLE_BBC(bX=RB$bX_RB,bXse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,cov_RB=RB$COV_RB,gcov=diag(p+1)*0,ldsc=rep(0,length(RB$by_RB)))
fit4=MRcML_BBC(bX=RB$bX_RB,bXse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,cov_RB=RB$COV_RB,gcov=diag(p+1)*0,ldsc=rep(0,length(RB$by_RB)))

#bX=RB$bX_RB;bXse=RB$bXse_RB;by=RB$by_RB;byse=RB$byse_RB;pv.thres=0.05;cov_RB=RB$COV_RB;gcov=diag(p+1)*0;ldsc=rep(0,length(RB$by_RB));max.iter=10;max.eps=1e-4;lambda=2;a=3;tau_upper=10;sampling.time=300;n_threads=4

BETA[i,,]=cbind(fit0$theta,fit1$theta,fit2$theta,fit3$theta,fit4$theta)
#print(i)
if(i%%10==0){
par(mfrow=c(2,2))
boxplot(BETA[1:i,1,])
lines(c(0:5),rep(theta0[1],6))
boxplot(BETA[1:i,2,])
lines(c(0:5),rep(theta0[2],6))
boxplot(BETA[1:i,3,])
lines(c(0:5),rep(theta0[3],6))
boxplot(BETA[1:i,4,])
lines(c(0:5),rep(theta0[4],6))
}
i=i+1
}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, " at ",i,"th iteration\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
i <<- i  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
}

