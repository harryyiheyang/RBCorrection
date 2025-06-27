#' Joint IV Selection with PLINK and Rao-Blackwell Correction for Multivariable MR
#'
#' Performs joint instrument selection using PLINK clumping based on a synthetic P-value
#' derived from the joint chi-square statistic across exposures. Applies Rao-Blackwellized
#' correction for selected SNP-exposure estimates.
#'
#' @param BETAMatrix A numeric matrix of SNP-exposure effect estimates (rows: SNPs, cols: exposures). Outcome should be the last column.
#' @param SEMatrix A numeric matrix of corresponding standard errors. Must match BETAMatrix in size and dimnames.Outcome should be the last column.
#' @param Rxy A (p+1 x p+1) covariance matrix of exposure effects. Outcome should be the last column and row.
#' @param bedfile Path to PLINK .bed/.bim/.fam file prefix.
#' @param p1 Clumping primary threshold (for selecting top SNPs).
#' @param p2 Clumping secondary threshold (for pruning LD neighbors).
#' @param r2 LD r-squared threshold.
#' @param kb Clumping distance (kb).
#' @param plink_path Full path to PLINK 1.9 binary (default = "./plink").
#' @param log_dir Optional path to save temporary clumping files. Uses a temp folder if not specified.
#' @param eta Standard deviation of noise for randomized Z (default = 1).
#' @param B Number of samples for Rao-Blackwell correction (default = 1000).
#'
#' @return A list with two matrices:
#' \describe{
#'   \item{\code{BETA_Sub}}{Rao-Blackwellized BETA matrix (rows: SNPs, cols: exposures).}
#'   \item{\code{SE_Sub}}{Total SE matrix (original SE + RB variance).}
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom data.table fwrite fread
#' @importFrom CppMatrix matrixVectorMultiply
#' @export
RB_MVJoint_Plink <- function(BETAMatrix, SEMatrix, Rxy, bedfile,
                             p1=5E-8, p2=5E-8, r2=1E-3, kb=1000,
                             plink_path = "./plink", log_dir = NULL,
                             eta = 1, B = 1000) {
  stopifnot(all(dim(BETAMatrix) == dim(SEMatrix)))
  rownames(BETAMatrix) <- rownames(BETAMatrix) %||% paste0("V", seq_len(nrow(BETAMatrix)))
  colnames(BETAMatrix) <- colnames(BETAMatrix) %||% paste0("V", seq_len(ncol(BETAMatrix)))
  rownames(SEMatrix) <- rownames(SEMatrix) %||% rownames(BETAMatrix)
  colnames(SEMatrix) <- colnames(SEMatrix) %||% colnames(BETAMatrix)
  stopifnot(all(rownames(BETAMatrix) == rownames(SEMatrix)))
  stopifnot(all(colnames(BETAMatrix) == colnames(SEMatrix)))

  SNPs <- rownames(BETAMatrix)
  p <- ncol(BETAMatrix) - 1
  m <- nrow(BETAMatrix)
  Thetaxx=solve(Rxy[1:p,1:p])

  # Step 1: randomized Z and joint chisq
  E = mvrnorm(n = m, mu = rep(0, p+1), Sigma = Rxy) * eta
  Z_base <- BETAMatrix / SEMatrix
  Z_rand <- Z_base + E
  chisq_stat <- apply(Z_rand[, 1:p, drop = FALSE], 1, function(z)
    sum(z * matrixVectorMultiply(Thetaxx, z))
  )
  P_joint <- pchisq(chisq_stat, df = p, lower.tail = FALSE)

  # Step 2: prepare PLINK clumping input
  gwas_df <- data.frame(SNP = SNPs, P = P_joint)
  if (is.null(log_dir)) {
    log_dir <- tempfile("plinkRBMVMR_")
    dir.create(log_dir)
    on.exit(unlink(log_dir, recursive = TRUE), add = TRUE)
  } else {
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  }
  gwas_path <- file.path(log_dir, "gwas.txt")
  fwrite(gwas_df, gwas_path, sep = "\t", quote = FALSE)

  # Step 3: switch directory to PLINK path
  old_dir <- getwd()
  plink_dir <- dirname(normalizePath(plink_path))
  setwd(plink_dir)
  on.exit(setwd(old_dir), add = TRUE)

  # Step 4: run PLINK clumping
  plink_cmd <- paste(
    "./plink",
    "--bfile", shQuote(bedfile),
    "--clump", shQuote(gwas_path),
    "--clump-p1", p1,
    "--clump-p2", p2,
    "--clump-r2", r2,
    "--clump-kb", kb,
    "--clump-field", "P",
    "--clump-snp-field", "SNP",
    "--out", file.path(log_dir, "clump")
  )
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Step 5: read selected SNPs
  clump_file <- file.path(log_dir, "clump.clumped")
  if (!file.exists(clump_file)) {
    warning("No SNPs passed PLINK clumping.")
    return(list(
      BETA_Sub = matrix(NA, nrow = 0, ncol = p+1),
      SE_Sub = matrix(NA, nrow = 0, ncol = p+1)
    ))
  }
  SNP_selected <- fread(clump_file)$SNP
  BETA_Select <- BETAMatrix[SNP_selected, , drop = FALSE]
  SE_Select <- SEMatrix[SNP_selected, , drop = FALSE]

  # Step 6: apply Rao-Blackwell correction
  cutoff <- qchisq(p1, df = p, lower.tail = FALSE)
  res <- RaoBlackwellCorrect_R(BETA_Select=BETA_Select, SE_Select=SE_Select,
                               Rxy=Rxy, eta=eta, cutoff=cutoff,B= B,onlyexposure=T)
  Beta_Sub <- res$BETA_RB
  SE_Sub <- res$SE_RB^2

  return(list(BETA_Sub = Beta_Sub, SE_Sub = SE_Sub))
}
