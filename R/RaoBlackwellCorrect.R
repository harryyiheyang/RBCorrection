#' Rao-Blackwell Correction for Selected IVs in MVMR
#'
#' Performs Rao-Blackwell correction via rejection sampling for each selected IV
#' based on randomized multivariate chi-square test.
#'
#' @param BETA_Select A numeric matrix of dimension m x p: effect estimates for selected IVs. Outcome should be the last column.
#' @param SE_Select A numeric matrix of same dimension: standard errors. Outcome should be the last column.
#' @param Rxy A (p+1) x (p+1) numeric matrix: covariance matrix of Z statistics. Outcome should be the last column and row.
#' @param eta Standard deviation of noise added to Z (default = 1).
#' @param pv.threshold P-value threshold used for selection.
#' @param B Number of samples used in rejection sampling (default = 1000).
#' @param onlyexposure A indicator of whether considering outcome when selecting IVs. Defaults to \code{TRUE}.
#' @param kappa_thres Numeric threshold controlling the removal of SNPs with unstable RBC covariance; if the condition number of the corrected covariance exceeds kappa_thres * kappa(Rxy) the SNP is discarded. Defaults to \code{10}.
#' @param warnings A indicator of whether printing warnings. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{BETA_RB}}{Corrected effect estimates (m x p).}
#'   \item{\code{SE_RB}}{Corrected standard errors (m x p).}
#' }
#'
#' @export
RaoBlackwellCorrect <- function(BETA_Select, SE_Select, Rxy, eta = 1, pv.threshold,
                                B = 1000, kappa_thres=10, onlyexposure=T,warnings=T) {
if(warnings){
cat("Please standardize data such that BETA = Zscore/sqrt n and SE = 1/sqrt n\n")
}
stopifnot(all(dim(BETA_Select) == dim(SE_Select)))
stopifnot(ncol(BETA_Select) == nrow(Rxy))
stopifnot(nrow(Rxy) == ncol(Rxy))

# Ensure matrix format with dimnames
p=ncol(BETA_Select)-1
if(onlyexposure==T){
cutoff=qchisq(pv.threshold,p,lower.tail=F)
}else{
cutoff=qchisq(pv.threshold,p+1,lower.tail=F)
}
BETA_Select <- as.matrix(BETA_Select)
SE_Select <- as.matrix(SE_Select)
Rxy <- as.matrix(Rxy)
Rxysqrt = matrixsqrt(Rxy)$w

if (is.null(rownames(BETA_Select))) rownames(BETA_Select) <- paste0("SNP", seq_len(nrow(BETA_Select)))
if (is.null(colnames(BETA_Select))) colnames(BETA_Select) <- paste0("Exp", seq_len(ncol(BETA_Select)))
rownames(SE_Select) <- rownames(BETA_Select)
colnames(SE_Select) <- colnames(BETA_Select)

res <- .Call(`_RBCorrection_RaoBlackwell`, beta_select=BETA_Select, se_select=SE_Select,
         Rxy=Rxy, Rxysqrt=Rxysqrt, eta=eta,cutoff=cutoff, B=B, onlyexposure=onlyexposure,n_threads=floor(parallel::detectCores() / 2))
res$SE_RB[res$CORRECTED_INDICES,]=sqrt(SE_Select[res$CORRECTED_INDICES,]^2*(1+1/eta^2))

for(i in 1:length(res$COV_RB)){
s=res$SE_RB[i,]
res$COV_RB[[i]]=t(t(Rxy)*s)*s-res$COV_RB[[i]]
s=sqrt(diag(res$COV_RB[[i]]))
s[is.na(s)]=0
res$SE_RB[i,]=s
}

condition_check = function(res, thres = 3) {
ind1 = which(!is.na(res$SE_RB[, p+1]) & res$SE_RB[, p+1] > 0)
ind2 = integer(0)
for (i in seq_along(res$COV_RB)) {
sevec   = SE_Select[i, ]
kapparxy = kappa(t(t(Rxy) * sevec) * sevec)
S        = res$COV_RB[[i]]
i1       = (diag(S) <= 0)
i2_flag  = kappa(S) > (thres * kapparxy)

if (sum(i1) + sum(i2_flag) == 0) {
ind2 = c(ind2, i)
}
}
ind = intersect(ind1, ind2)
return(ind)
}


ind=condition_check(res,thres=kappa_thres)
res$BETA_RB=res$BETA_RB[ind,]
res$SE_RB=res$SE_RB[ind,]
res$COV_RB=res$COV_RB[ind]

rownames(res$BETA_RB) <- rownames(BETA_Select[ind,])
colnames(res$BETA_RB) <- colnames(BETA_Select)
rownames(res$SE_RB) <- rownames(BETA_Select[ind,])
colnames(res$SE_RB) <- colnames(BETA_Select)

bX_RB=res$BETA_RB[,1:p]
bXse_RB=res$SE_RB[,1:p]
by_RB=res$BETA_RB[,1+p]
byse_RB=res$SE_RB[,1+p]
if(p==1){
bX_RB=as.vector(bX_RB)
bXse_RB=as.vector(bXse_RB)
}
return(list(bX_RB=bX_RB,bXse_RB=bXse_RB,by_RB=by_RB,byse_RB=byse_RB,COV_RB=res$COV_RB))
}
