library(MRBEE)
library(MASS)
library(devtools)
library(MRAPSS)
library(dplyr)
document()
m=1000
n=1e4
Rxy=(matrix(0.5,2,2)+diag(2)/2)
#Rxy=diag(2)
BETA=matrix(0,1000,5)
source("example/otherfunction.R")

i=1
while(i<=1000){
indicator=F
tryCatch({
betax=rnorm(n=m,mean=0,sd=1)
betax[sample(m,0.8*m)]=0
betax=betax/sqrt(sum(betax^2))*sqrt(0.05)
betay=betax*0.2
E=mvrnorm(n=m,mu=rep(0,2),Sigma=Rxy)
E=E%*%diag(sqrt(c(1/n,1/n)))
bx=betax+E[,1]
by=betay+E[,2]
bxse=byse=rep(sqrt(1/n),m)
pv=pchisq((bx/bxse)^2,1,lower.tail=F)
indselect=which(pv<0.05)
MRData=data.frame(SNP=paste0("V",1:m),A1=rep("T",m),A2=rep("G",m),b.exp=bx,b.out=by,se.exp=bxse,se.out=byse)%>%mutate(pval.exp=pchisq(b.exp^2/se.exp^2,1,lower.tail=F),pval.out=pchisq(b.out^2/se.out^2,1,lower.tail=F),L2=1,Threshold=0.05)
fitAPSS=MRAPSS(MRData,exposure = "X",outcome = "Y",C = Rxy,Omega =  Rxy*0 ,Cor.SelectionBias = T,tol=1e-5)
fit0=MRBEE.IMRP.UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy,pv.thres=0)
fit1=MRBEE.IMRP.UV(by=by[indselect],bx=bx[indselect],byse=byse[indselect],bxse=bxse[indselect],Rxy=Rxy,pv.thres=0)
RB=RB_UV_Analytic(gamma=bx,sigma=bxse,cutoff=sqrt(qchisq(0.05,1,lower.tail=F)),eta=1)
fit2=MRBEE.IMRP.UV(by=by[RB$IVselect],byse=byse[RB$IVselect],bx=RB$BETA_RB,bxse=RB$SE_RB,Rxy=Rxy,pv.thres=0)
RB=RB_UV_Simulation_Joint(GWAS_X=data.frame(SNP=paste0("V",1:m),BETA=bx,SE=bxse),
                          GWAS_Y=data.frame(SNP=paste0("V",1:m),BETA=by,SE=byse),
                          Rxy=Rxy,P_threshold=0.05,eta=1.5,B=1500)
fit3=MRBEE_UV_Winner(bx=RB$GWAS_X_RB$BETA_RB,bxse=RB$GWAS_X_RB$SE_RB,by=RB$GWAS_Y_RB$BETA_RB,byse=RB$GWAS_Y_RB$SE_RB,Rxy=Rxy,pv.thres=0,S_RB=RB$S_RB,phi=1)
BETA[i,]=c(fitAPSS$beta,fit0$theta,fit1$theta,fit2$theta,fit3$theta)
i=i+1
print(i)
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
boxplot(BETA[1:277,],ylim=c(-1,1))
lines(c(0:10),rep(0.2,11))
