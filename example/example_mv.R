library(MRBEE)
library(MASS)
library(devtools)
library(MRAPSS)
library(dplyr)
document()
m=1000
n=3e4
p=4
Rbb=matrix(0.5,4,4)+diag(4)*0.5
Rxy=matrix(0.9,5,5)+diag(5)*0.1
Rxy=diag(5)
Thetaxx=solve(Rxy[1:p,1:p])
BETA=array(0,c(1000,4,3))
theta0=c(0,0.2,0,0)
iv.thes=5e-3
h2=0.05
i=1

while(i<=1000){
indicator=F
tryCatch({
betaX=mvrnorm(n=m,mu=rep(0,p),Rbb)
for(ii in 1:p){
betax=betaX[,ii]
#ind=sample(m,0.5*m)
#betax[ind]=0
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
RB=RB_MV_Joint(BETAMatrix=BETAMatrix,SEMatrix=SEMatrix,
                          Rxy=Rxy,P_threshold=iv.thes,eta=1,B=1000)
print(c(length(indselect),dim(RB$BETA_RB)[1]))
fit2=MRBEE_Winner(bX=RB$BETA_RB[,1:p],bXse=RB$SE_RB[,1:p],by=RB$BETA_RB[,p+1],byse=RB$SE_RB[,p+1],Rxy=Rxy,pv.thres=0.05,S_RB=RB$COV_RB,phi=1)
BETA[i,,]=cbind(fit0$theta,fit1$theta,fit2$theta)
#print(i)
if(i%%10==0){
par(mfrow=c(2,2))
boxplot(BETA[1:i,1,])
lines(c(0:4),rep(theta0[1],5))
boxplot(BETA[1:i,2,])
lines(c(0:4),rep(theta0[2],5))
boxplot(BETA[1:i,3,])
lines(c(0:4),rep(theta0[3],5))
boxplot(BETA[1:i,4,])
lines(c(0:4),rep(theta0[4],5))
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

