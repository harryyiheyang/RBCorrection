#' Joint Rao-Blackwellized Correction for Univariate MR (Exposure + Outcome)
#'
#' @param GWAS_X Data.frame with columns "SNP", "BETA", "SE" for the exposure.
#' @param GWAS_Y Data.frame with columns "SNP", "BETA", "SE" for the outcome, rows must match GWAS_X.
#' @param Rxy 2x2 correlation/covariance matrix for estimation errors (rows/columns: exposure and outcome).
#' @param bedfile Path to PLINK .bed/.bim/.fam prefix.
#' @param p1,p2,r2,kb PLINK clumping parameters.
#' @param plink_path Path to PLINK binary.
#' @param log_dir Temporary directory to store intermediate files.
#' @param eta Noise sd for randomized selection Z.
#' @param B Number of Monte Carlo samples in rejection sampling.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats qchisq rnorm dnorm pnorm
#' @return A list with two data.frames: GWAS_X_RB and GWAS_Y_RB
#' @export
RB_UV_Plink_Joint <- function(GWAS_X, GWAS_Y, Rxy,
                              bedfile,
                              p1 = 5E-8, p2 = 5E-8, r2 = 1E-3, kb = 1000,
                              plink_path = "./plink",
                              log_dir = NULL, eta = 1, B = 1000) {

  stopifnot(all(c("SNP", "BETA", "SE") %in% names(GWAS_X)))
  stopifnot(all(c("SNP", "BETA", "SE") %in% names(GWAS_Y)))
  stopifnot(nrow(GWAS_X) == nrow(GWAS_Y))
  stopifnot(all(GWAS_X$SNP == GWAS_Y$SNP))
  stopifnot(all(dim(Rxy) == c(2, 2)))

  # Step 1: Create log directory
  if (is.null(log_dir)) {
    log_dir <- tempfile("plinkRB_")
    dir.create(log_dir)
    on.exit(unlink(log_dir, recursive = TRUE), add = TRUE)
  } else {
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Step 2: Generate randomized Z for X
  n <- nrow(GWAS_X)
  cutoff <- sqrt(qchisq(p1, df = 1, lower.tail = FALSE))
  z_base <- GWAS_X$BETA / GWAS_X$SE
  e = mvrnorm(n,rep(0,2),Rxy)*eta
  z_rand <- z_base + e[,1]

  # Pre-filter for PLINK input
  sel <- abs(z_rand) > (cutoff / 1.8)
  if (sum(sel) == 0) {
    warning("No SNPs passed randomized pre-selection for clumping.")
    return(list(
      GWAS_X_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric()),
      GWAS_Y_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric())
    ))
  }

  clump_input <- GWAS_X[sel, ]
  chisq_stat <- (clump_input$BETA / clump_input$SE)^2
  clump_input$P <- pchisq(chisq_stat, df = 1, lower.tail = FALSE)

  gwas_path <- file.path(log_dir, "gwas.txt")
  data.table::fwrite(clump_input[, c("SNP", "P")], gwas_path,
                     sep = "\t", row.names = FALSE, quote = FALSE)

  # Step 3: Run PLINK clumping
  old_dir <- getwd()
  plink_dir <- dirname(normalizePath(plink_path))
  setwd(plink_dir)
  on.exit(setwd(old_dir), add = TRUE)

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

  # Step 4: Load clumped SNPs
  clump_file <- file.path(log_dir, "clump.clumped")
  if (!file.exists(clump_file)) {
    warning("PLINK did not return clumped SNPs.")
    return(list(
      GWAS_X_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric()),
      GWAS_Y_RB = data.frame(SNP = character(), BETA_RB = numeric(), SE_RB = numeric())
    ))
  }

  SNPs <- data.table::fread(clump_file)$SNP
  sel <- GWAS_X$SNP %in% SNPs

  # Step 5: Perform Joint RB correction
  gamma <- GWAS_X$BETA[sel]
  Gamma <- GWAS_Y$BETA[sel]
  sigma_X <- GWAS_X$SE[sel]
  sigma_Y <- GWAS_Y$SE[sel]
  BETA_Select=cbind(gamma,Gamma)
  SE_Select=cbind(sigma_X,sigma_Y)
  RB <- RaoBlackwellCorrect_R(BETA_Select=BETA_Select, SE_Select=SE_Select, Rxy=Rxy,
                              eta=eta, cutoff=cutoff^2,B= B,onlyexposure=T)

  out_X <- data.frame(SNP = GWAS_X$SNP[sel],
                      BETA_RB = RB$BETA_RB[,1],
                      SE_RB = RB$SE_RB[,1])

  out_Y <- data.frame(SNP = GWAS_Y$SNP[sel],
                      BETA_RB = RB$BETA_RB[,2],
                      SE_RB = RB$SE_RB[,2])

  return(list(GWAS_X_RB = out_X, GWAS_Y_RB = out_Y))
}
