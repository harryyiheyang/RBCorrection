covfreq=function(theta,theta.se,theta0){
g=ifelse(abs(theta/1.96-theta0/1.96)<=theta.se,1,0)
return(g)
}

rejfreq=function(theta,theta.se){
g=ifelse(abs(theta/1.96)<=theta.se,0,1)
return(g)
}
library(MRBEE)
library(MASS)
library(devtools)
library(MRAPSS)
library(dplyr)
document()
devtools::load_all("C:/Users/yxy1234/Downloads/MRcare-main/")
m=500
n=10e4
Rxy=(matrix(0.5,2,2)+diag(2)/2)
#Rxy=diag(2)
BETA=SE=COV=REJ=matrix(0,300,8)
colnames(BETA)=c("APSS","CARE","BEE-Causal","BEE-Naive","BEE-RBC","BEE-RBS","RAPS-RBS","cML-RBS")
theta0=0.5
IV.theshold=5e-8
pleio.ratio=0
eta=0.5

i=1
while(i<=300){
indicator=F
tryCatch({
betax=rnorm(n=m,mean=0,sd=1)
betax[sample(m,50)]=0
betax=betax/sqrt(sum(betax^2))*sqrt(0.05)
betay=betax*theta0
E=mvrnorm(n=m,mu=rep(0,2),Sigma=Rxy)
E=E%*%diag(sqrt(c(1/n,1/n)))
bx=betax+E[,1]
by=betay+E[,2]+ifelse(pleio.ratio>0,rbinom(m,1,pleio.ratio)*rnorm(m,0,sqrt(0.001/pleio.ratio)),0)
bxse=byse=rep(sqrt(1/n),m)
print((mean(bx^2)-bxse[1]^2)/mean(bx^2))
pv=pchisq((bx/bxse)^2,1,lower.tail=F)
indselect=which(pv<IV.theshold)
MRData=data.frame(SNP=paste0("V",1:length(indselect)),A1=rep("T",length(indselect)),A2=rep("G",length(indselect)),b.exp=bx[indselect],b.out=by[indselect],se.exp=bxse[indselect],se.out=byse[indselect])%>%mutate(pval.exp=pchisq(b.exp^2/se.exp^2,1,lower.tail=F),pval.out=pchisq(b.out^2/se.out^2,1,lower.tail=F),L2=1,Threshold=IV.theshold)
fitAPSS=MRAPSS(MRData,exposure = "X",outcome = "Y",C = Rxy,Omega =  Rxy*0 ,Cor.SelectionBias = T,tol=1e-5)
fitCARE=CARE2_boot(gamma.exp_sel=bx[indselect],gamma.out_sel=by[indselect],se.exp_sel=bxse[indselect],se.out_sel=byse[indselect],nx=n,ny=n,nrep=100,biascorrect="rerand",algorithm="CD",pthr=IV.theshold)
fit0=MRBEE.IMRP.UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy,pv.thres=0.05)
fit1=MRBEE.IMRP.UV(by=by[indselect],bx=bx[indselect],byse=byse[indselect],bxse=bxse[indselect],Rxy=Rxy,pv.thres=0.05)

RB=RB_UV_Analytic(gamma=bx,sigma=bxse,cutoff=sqrt(qchisq(IV.theshold,1,lower.tail=F)),eta=eta)
fit2=MRBEE.IMRP.UV(by=by[RB$IVselect],byse=byse[RB$IVselect],bx=RB$BETA_RB,bxse=RB$SE_RB,Rxy=Rxy,pv.thres=0.05)

pv=pchisq((bx/bxse)^2+rnorm(m,0,eta)^2,1,lower.tail=F)
indselect=which(pv<IV.theshold)
RB=RaoBlackwellCorrect(BETA_Select=cbind(bx[indselect],by[indselect]),SE_Select=cbind(bxse[indselect],byse[indselect]),Rxy=Rxy,eta=eta,pv.threshold=IV.theshold,B=1000,warnings=F)
fit3=MRBEE_UV_BBC(bx=RB$bX_RB,bxse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,pv.thres=0.05,cov_RB=RB$COV_RB,gcov=diag(2)*0,sampling.strategy="bootstrap")
fit4=MRRAPS_BBC(bx=RB$bX_RB,bxse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,cov_RB=RB$COV_RB,gcov=diag(2)*0,sampling.strategy="bootstrap")
fit5=MRcML_UV_BBC(bx=RB$bX_RB,bxse=RB$bXse_RB,by=RB$by_RB,byse=RB$byse_RB,cov_RB=RB$COV_RB,gcov=diag(2)*0,sampling.strategy="bootstrap")

BETA[i,]=c(fitAPSS$beta,fitCARE$res$MA_BIC["tilde_theta"],fit0$theta,fit1$theta,fit2$theta,fit3$theta,fit4$theta,fit5$theta)
SE[i,]=sqrt(c(fitAPSS$beta.se^2,fitCARE$res$MA_BIC["Efron_se"]^2,fit0$vartheta,fit1$vartheta,fit2$vartheta,fit3$theta.var,fit4$theta.var,fit5$theta.var))
COV[i,]=covfreq(BETA[i,],SE[i,],theta0)
REJ[i,]=rejfreq(BETA[i,],SE[i,])
i=i+1
if(i%%10==0){
boxplot(BETA[1:i,],ylim=c(0.4,0.6))
lines(c(0:8),rep(theta0,9))
print(colMeans(COV[1:i,]))
}
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
boxplot(BETA)
lines(c(0:10),rep(0.2,11))
