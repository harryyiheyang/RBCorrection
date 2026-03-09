#!/usr/bin/env Rscript

## =========================================================
## MV simulation (parallel)
## = your old MV DGP + UV-style future_lapply wrapper
## + mixture (small/med/large) sparse effects
## + sparse heritability targets (ABSOLUTE), while KEEPING
##   original background heritability (Sigma_polygenic uses h2x/h2y as before)
## =========================================================

suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(MRcML)
  library(MVMRcML)
  library(MRBEE)
  library(MRcare)
  library(MASS)
  library(mr.raps)
  library(CppMatrix)
  library(Rcpp)
  library(RBCorrection)
  library(future.apply)
})

## RBCorrection analytic helper (defines rb_onlyexposure_analytic_batch)
RB_HELPER <- "/projects/standard/panwei/lin00374/RBCorrection/R/RB_analytic.R"
source(RB_HELPER)
## -------------------------
## Mixture MV generator
## -------------------------
rmvnorm_mixture <- function(n, Sigma_base,
                            w = c(0.85, 0.14, 0.01),
                            sd_rel = c(1, 4, 12)) {
  stopifnot(n >= 0, length(w) == length(sd_rel))
  if (n == 0) return(matrix(numeric(0), nrow = 0, ncol = ncol(Sigma_base)))

  w <- w / sum(w)
  comp <- sample.int(length(w), size = n, replace = TRUE, prob = w)
  scale <- sd_rel[comp]

  U <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Sigma_base)), Sigma = Sigma_base)
  sweep(U, 1, scale, "*")
}

rescale_matrix_col_h2 <- function(G, h2_per_col) {
  # G: n x K
  # h2_per_col: either scalar or length-K vector
  K <- ncol(G)
  if (length(h2_per_col) == 1) h2_per_col <- rep(h2_per_col, K)
  stopifnot(length(h2_per_col) == K)

  ss <- colSums(G^2)
  scale <- rep(1, K)
  ok <- is.finite(ss) & ss > 0 & is.finite(h2_per_col) & h2_per_col >= 0
  scale[ok] <- sqrt(h2_per_col[ok] / ss[ok])

  sweep(G, 2, scale, "*")
}

rescale_vector_h2 <- function(x, h2_target) {
  ss <- sum(x^2)
  if (!is.finite(ss) || ss <= 0 || !is.finite(h2_target) || h2_target < 0) return(x)
  x * sqrt(h2_target / ss)
}

## -------------------------
## Arguments
## (same 9 args as old MV, plus OPTIONAL:
##  10) h2_gamma_perX (absolute, per exposure)
##  11) h2_xi         (absolute, direct pleiotropy in Y)
## )
## -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Need at least 9 args: indx1 Nin PropInvalidIn h2x h2y rhoxy rhoos rhoos_x job.id [h2_gamma_perX] [h2_xi]")
}

indx1         <- as.numeric(args[[1]])
Nin           <- as.character(args[[2]])
PropInvalidIn <- as.character(args[[3]])
h2x           <- as.numeric(args[[4]])
h2y           <- as.numeric(args[[5]])
rhoxy         <- as.numeric(args[[6]])
rhoos         <- as.numeric(args[[7]])
rhoos_x       <- as.numeric(args[[8]])
job.id        <- as.numeric(args[[9]])

## OPTIONAL targets (absolute)
## If you don’t pass them, these defaults will be used.
h2_gamma_perX <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.10
h2_xi_target  <- if (length(args) >= 11) as.numeric(args[[11]]) else -9

if (!is.finite(h2_gamma_perX) || h2_gamma_perX < 0) stop("h2_gamma_perX must be >= 0")

indx3 <- match(Nin, c("N1","N2","N3","N4","N5"))
if (is.na(indx3)) stop("Nin must be one of N1..N5")

indx4 <- match(PropInvalidIn, c("Prop1","Prop2","Prop3"))
if (is.na(indx4)) stop("PropInvalidIn must be one of Prop1..Prop3")

## -------------------------
## Main settings
## -------------------------
save_datdir <- "/projects/standard/panwei/lin00374/RBCorrection/Simulation_mixture/MV/"
if (!dir.exists(save_datdir)) dir.create(save_datdir, recursive = TRUE)

eta <- 0.5

thetavec <- rbind(
  c(0,   0.2,    0),
  c(0.1, 0.05, 0)
)
Nvec <- c(5e4, 8e4, 1e5, 2.5e5, 5e5)
prop_invalid_vec <- c(0, 0.1, 0.3)

theta <- thetavec[indx1, ]
N <- Nvec[indx3]
prop_invalid <- prop_invalid_vec[indx4]

pthr  <- 5e-8
pthr2 <- 5e-5
NxNy_ratio <- 1

M <- 5e4
K <- 3

pi1 <- 0.01 * (1 - prop_invalid)
pi2 <- 0.01 * prop_invalid

nx <- N
ny <- N / NxNy_ratio

## -------------------------
## Covariance structures (same as your old MV code)
## -------------------------
rhoxx <- rbb <- rhoxy
Rbb <- matrix(rbb, nrow = K, ncol = K); diag(Rbb) <- 1

sigma2x      <- 1e-4
sigma2x_td   <- 1e-4
sigma2y_td   <- 1e-4

Sigma_b_valid <- sigma2x    * Rbb
Sigma_b_inval <- sigma2x_td * Rbb

## polygenic background = ORIGINAL (uses h2x/h2y directly, no subtraction)
Sigma_polygenic <- matrix(0, nrow = K + 1, ncol = K + 1)
Sigma_polygenic[1:K, 1:K] <- rhoxx * h2x
diag(Sigma_polygenic)[1:K] <- h2x
Sigma_polygenic[K+1, K+1] <- h2y
Sigma_polygenic[1:K, K+1] <- rhoxy * sqrt(h2x) * sqrt(h2y)
Sigma_polygenic[K+1, 1:K] <- rhoxy * sqrt(h2x) * sqrt(h2y)
Sigma_polygenic <- Sigma_polygenic / (7e6)

## sampling covariance (same as old)
Sigma_delta <- matrix(0, nrow = K + 1, ncol = K + 1)
Sigma_delta[1:K, 1:K] <- rhoos_x * (sqrt(1/nx)^2)
diag(Sigma_delta)[1:K] <- sqrt(1/nx)^2
Sigma_delta[1:K, K+1] <- rhoos * sqrt(1/nx) * sqrt(1/ny)
Sigma_delta[K+1, 1:K] <- rhoos * sqrt(1/nx) * sqrt(1/ny)
Sigma_delta[K+1, K+1] <- sqrt(1/ny)^2

Thetaxx <- solve(cov2cor(Sigma_delta[1:K, 1:K]))

Vx <- Sigma_delta[1:K, 1:K] + Sigma_polygenic[1:K, 1:K]
Vx_inv <- solve(Vx)


cat(
  paste(
    "N", N,
    "pthr", pthr,
    "theta", paste(theta, collapse = ","),
    "prop_invalid", prop_invalid,
    "rhoxy", rhoxy,
    "rhoos", rhoos,
    "rhoos_x", rhoos_x,
    "| background h2x", h2x, "h2y", h2y,
    "| sparse targets: h2_gamma_perX", h2_gamma_perX, "h2_xi", h2_xi_target
  ),
  "\n"
)

## -------------------------
## Simulation index set
## -------------------------
set.ind <- job.id
simulation.ind.set <- ((set.ind - 1) * 100 + 1):(set.ind * 100)

## -------------------------
## One replicate
## -------------------------
run_one_sim <- function(sim.ind,
                        M, K, pi1, pi2, theta,
                        Sigma_b_valid, Sigma_b_inval, sigma2y_td,
                        Sigma_polygenic, Sigma_delta,
                        nx, ny,
                        pthr, pthr2,
                        Thetaxx, Vx, Vx_inv,
                        eta,
                        h2_gamma_perX, h2_xi_target,
                        RB_HELPER) {

  ## RB helper must exist in each worker (important for multisession)
  if (!exists("rb_onlyexposure_analytic_batch", mode = "function")) {
    source(RB_HELPER)
  }

  cat("simulation:", sim.ind, "\n")

  numIV <- 0
  numIV2 <- 0
  tmpj <- 0

  setting <- NULL

  ## placeholders that must exist after loop
  ind_filter <- integer(0)
  betahat_x <- betahat_y <- NULL
  se_x <- se_y <- NULL

  while (numIV < 3 || numIV2 < 4) {

    set.seed(sim.ind + 10000 * tmpj)

    ind1 <- sample.int(M, round(M * pi1))
    causalsnps <- ind1
    ind2 <- sample(setdiff(seq_len(M), causalsnps), round(M * pi2))
    causalsnps <- c(causalsnps, ind2)

    gamma <- matrix(0, nrow = M, ncol = K)
    xi <- rep(0, M)

    ## mixture sparse effects
    if (length(ind1) > 0) {
      gamma[ind1, ] <- rmvnorm_mixture(n = length(ind1), Sigma_base = Sigma_b_valid)
    }
    if (length(ind2) > 0) {
      gamma[ind2, ] <- rmvnorm_mixture(n = length(ind2), Sigma_base = Sigma_b_inval)

      ## xi initial draw (will rescale to target)
      xi[ind2] <- rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
    }

    ## rescale sparse components to hit target h2 (absolute)

    if (length(causalsnps) > 0) {
      gamma[causalsnps, ] <- rescale_matrix_col_h2(
        gamma[causalsnps, , drop = FALSE],
        h2_per_col = h2_gamma_perX
      )
    }
    if (length(ind2) > 0 & h2_xi_target > 0) {
      xi[ind2] <- rescale_vector_h2(xi[ind2], h2_xi_target)
    }

    ## polygenic background (unchanged)
    poly_effect <- MASS::mvrnorm(n = M, mu = rep(0, K + 1), Sigma = Sigma_polygenic)

    ## true marginal effects
    betax <- gamma + poly_effect[, 1:K, drop = FALSE]
    betay <- as.vector(gamma %*% theta) + xi + poly_effect[, K + 1]

    ## sampling error
    delta <- MASS::mvrnorm(n = M, mu = rep(0, K + 1), Sigma = Sigma_delta)

    betahat_x <- betax + delta[, 1:K, drop = FALSE]
    betahat_y <- betay + delta[, K + 1]

    se_x <- matrix(sqrt(1 / nx), nrow = M, ncol = K)
    se_y <- rep(sqrt(1 / ny), M)

    ## MV selection (same as your old MV)
    zx_mat <- betahat_x / se_x
    chi2_x <- rowSums((zx_mat %*% Thetaxx) * zx_mat)
    pv_x <- pchisq(chi2_x, df = K, lower.tail = FALSE)

    ind_filter   <- which(pv_x < pthr)
    ind_filter_2 <- which(pv_x < pthr2)

    numIV  <- length(ind_filter)
    numIV2 <- length(ind_filter_2)

    ## diagnostics
    hertx <- colSums(betax^2)
    herty <- sum(betay^2)
    hertx.select <- if (numIV > 0) colSums(betax[ind_filter, , drop = FALSE]^2) else rep(NA_real_, K)

    Fstat <- if (numIV > 0) {
      colMeans(betahat_x[ind_filter, , drop = FALSE]^2 / se_x[ind_filter, , drop = FALSE]^2) - 1
    } else rep(NA_real_, K)

    varX <- if (numIV > 0) colSums(betax[ind_filter, , drop = FALSE]^2) else rep(NA_real_, K)
    varY <- if (numIV > 0) sum(betay[ind_filter]^2) else NA_real_

    nind1 <- sum(ind_filter %in% ind1)
    nind2 <- sum(ind_filter %in% ind2)

    sim.setting <- c(
      nIV = numIV, nind1 = nind1, nind2 = nind2,
      hertx, hertx.select, herty,
      Fstat,
      varX, varY
    )
    names(sim.setting) <- c(
      "nIV","nind1","nind2",
      paste0("hertX", 1:K), paste0("hertX.select", 1:K), "hertY",
      paste0("F", 1:K),
      paste0("varX", 1:K), "varY"
    )
    setting <- sim.setting

    tmpj <- tmpj + 1
    cat("numIV:", numIV, "numIV2:", numIV2, "\n")
  }

  ## -----------------------------
  ## RBCorrect block (your old MV)
  ## -----------------------------
  BEE_RB.est <- cML_RB.est <- MRcML_RB.est <- rep(NA_real_, K)
  BEE_RB.se  <- cML_RB.se  <- MRcML_RB.se  <- rep(NA_real_, K)

  run.time <- c(BEE_RB = NA_real_, cML_RB = NA_real_, MRcML_RB = NA_real_,
                BEE = NA_real_, MRcML = NA_real_, IVW = NA_real_, Egger = NA_real_, Median = NA_real_)

  Z_beta <- MASS::mvrnorm(n = M, mu = rep(0, K), Sigma = (eta^2) * Vx)
  b_tilde <- betahat_x + Z_beta

  Sstat <- rowSums((b_tilde %*% Vx_inv) * b_tilde)
  pv_rb <- pchisq(Sstat, df = K, lower.tail = FALSE)
  indselect <- which(pv_rb < pthr)

  if (length(indselect) > 2) {

    RB <- tryCatch(
      rb_onlyexposure_analytic_batch(
        beta_select = cbind(betahat_x[indselect, , drop = FALSE], betahat_y[indselect]),
        SE_Select   = cbind(se_x[indselect, , drop = FALSE], se_y[indselect]),
        gcov = Sigma_polygenic,
        ldsc = rep(1, length(indselect)),
        Rxy  = cov2cor(Sigma_delta),
        eta  = eta,
        pv.threshold = pthr
      ),
      error = function(e) { cat("ERROR RB:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(RB)) {

      t0 <- proc.time()[3]
      BEE_RB <- tryCatch(
        MRBEE_BBC(
          bX = RB$bX_RB, bXse = RB$bXse_RB,
          by = RB$by_RB, byse = RB$byse_RB,
          pv.thres = 0.05,
          cov_RB = RB$COV_RB,
          sampling.strategy = "bootstrap",
          n_threads = 1
        ),
        error = function(e) { cat("ERROR BEE_RB:", conditionMessage(e), "\n"); NULL }
      )
      run.time["BEE_RB"] <- proc.time()[3] - t0
      BEE_RB.est <- if (is.null(BEE_RB)) rep(NA_real_, K) else as.numeric(BEE_RB$theta)
      BEE_RB.se  <- if (is.null(BEE_RB)) rep(NA_real_, K) else as.numeric(BEE_RB$theta.se)

      t0 <- proc.time()[3]
      cML_RB <- tryCatch(
        MRcML_BBC(
          bX = RB$bX_RB, bXse = RB$bXse_RB,
          by = RB$by_RB, byse = RB$byse_RB,
          cov_RB = RB$COV_RB,
          sampling.strategy = "bootstrap",
          n_threads = 1
        ),
        error = function(e) { cat("ERROR cML_RB:", conditionMessage(e), "\n"); NULL }
      )
      run.time["cML_RB"] <- proc.time()[3] - t0
      cML_RB.est <- if (is.null(cML_RB)) rep(NA_real_, K) else as.numeric(cML_RB$theta)
      cML_RB.se  <- if (is.null(cML_RB)) rep(NA_real_, K) else as.numeric(cML_RB$theta.se)

      Glist <- tryCatch(lapply(RB$COV_RB, solve), error = function(e) NULL)
      if (!is.null(Glist) && !is.null(cML_RB)) {
        K_center <- sum(cML_RB$gamma != 0)
        K_lo <- pmax(0, K_center - 5)
        K_hi <- pmin(K_center + 5, length(RB$bX_RB) - 2)
        if (K_lo <= K_hi) {
          t0 <- proc.time()[3]
          MRcML_RB <- tryCatch(
            MVMRcML::MVmr_cML_DP(
              b_exp = as.matrix(RB$bX_RB),
              b_out = as.matrix(RB$by_RB),
              se_bx = as.matrix(RB$bXse_RB),
              Sig_inv_l = Glist,
              n = min(nx, ny),
              K_vec = K_lo:K_hi,
              maxit = 1000,
              thres = 1e-4
            ),
            error = function(e) { cat("ERROR MRcML_RB:", conditionMessage(e), "\n"); NULL }
          )
          run.time["MRcML_RB"] <- proc.time()[3] - t0
          MRcML_RB.est <- if (is.null(MRcML_RB)) rep(NA_real_, K) else as.numeric(MRcML_RB$BIC_DP_theta)
          MRcML_RB.se  <- if (is.null(MRcML_RB)) rep(NA_real_, K) else as.numeric(MRcML_RB$BIC_DP_se)
        }
      }
    }
  }

  ## -----------------------------
  ## Other methods (same as old MV)
  ## -----------------------------
  BEE.est <- MRcML.est <- IVW.est <- Egger.est <- Median.est <- rep(NA_real_, K)
  BEE.se  <- MRcML.se  <- IVW.se  <- Egger.se  <- Median.se  <- rep(NA_real_, K)

  if (length(ind_filter) > 3) {

    b_exp  <- betahat_x[ind_filter, , drop = FALSE]
    b_out  <- betahat_y[ind_filter]
    se_exp <- se_x[ind_filter, , drop = FALSE]
    se_out <- se_y[ind_filter]

    t0 <- proc.time()[3]
    BEE.result <- tryCatch(
      MRBEE::MRBEE.IMRP(
        by = b_out,
        bX = as.matrix(b_exp),
        byse = se_out,
        bXse = as.matrix(se_exp),
        Rxy = cov2cor(Sigma_delta)
      ),
      error = function(e) { cat("ERROR BEE:", conditionMessage(e), "\n"); NULL }
    )
    run.time["BEE"] <- proc.time()[3] - t0
    BEE.est <- if (is.null(BEE.result)) rep(NA_real_, K) else as.numeric(BEE.result$theta)
    BEE.se  <- if (is.null(BEE.result)) rep(NA_real_, K) else as.numeric(sqrt(diag(BEE.result$covtheta)))

    t0 <- proc.time()[3]
    Sig_inv_l <- MVMRcML::invcov_mvmr(se_bx = as.matrix(se_exp), se_by = se_out, rho_mat = cov2cor(Sigma_delta))

    res0 <- tryCatch(
      MVMRcML::MVmr_cML(
        b_exp = as.matrix(b_exp),
        b_out = as.matrix(b_out),
        se_bx = as.matrix(se_exp),
        Sig_inv_l = Sig_inv_l,
        n = min(nx, ny),
        maxit = 1000,
        thres = 1e-4
      ),
      error = function(e) { cat("ERROR MVmr_cML:", conditionMessage(e), "\n"); NULL }
    )

    K_vec <- 0:min(length(ind2), length(b_out) - K - 1)
    if (!is.null(res0) && !is.null(res0$BIC_invalid)) {
      K_vec <- (pmax(0, length(res0$BIC_invalid) - 5)):(pmin(length(res0$BIC_invalid) + 5, length(b_out) - K - 1))
    }

    res_MRcML <- tryCatch(
      MVMRcML::MVmr_cML_DP(
        b_exp = as.matrix(b_exp),
        b_out = as.matrix(b_out),
        se_bx = as.matrix(se_exp),
        Sig_inv_l = Sig_inv_l,
        n = min(nx, ny),
        K_vec = K_vec,
        maxit = 1000,
        thres = 1e-4
      ),
      error = function(e) { cat("ERROR MVmr_cML_DP:", conditionMessage(e), "\n"); NULL }
    )
    run.time["MRcML"] <- proc.time()[3] - t0
    MRcML.est <- if (is.null(res_MRcML)) rep(NA_real_, K) else as.numeric(res_MRcML$BIC_DP_theta)
    MRcML.se  <- if (is.null(res_MRcML)) rep(NA_real_, K) else as.numeric(res_MRcML$BIC_DP_se)

    mvobj <- MendelianRandomization::mr_mvinput(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)

    t0 <- proc.time()[3]
    set.seed(1)
    MR_IVW <- tryCatch(MendelianRandomization::mr_mvivw(mvobj), error = function(e) NULL)
    run.time["IVW"] <- proc.time()[3] - t0
    IVW.est <- if (is.null(MR_IVW)) rep(NA_real_, K) else as.numeric(MR_IVW@Estimate)
    IVW.se  <- if (is.null(MR_IVW)) rep(NA_real_, K) else as.numeric(MR_IVW@StdError)

    t0 <- proc.time()[3]
    set.seed(1)
    MR_EGGER <- tryCatch(MendelianRandomization::mr_mvegger(mvobj), error = function(e) NULL)
    run.time["Egger"] <- proc.time()[3] - t0
    Egger.est <- if (is.null(MR_EGGER)) rep(NA_real_, K) else as.numeric(MR_EGGER@Estimate)
    Egger.se  <- if (is.null(MR_EGGER)) rep(NA_real_, K) else as.numeric(MR_EGGER@StdError.Est)

    t0 <- proc.time()[3]
    set.seed(1)
    MR_MED <- tryCatch(MendelianRandomization::mr_mvmedian(mvobj), error = function(e) NULL)
    run.time["Median"] <- proc.time()[3] - t0
    Median.est <- if (is.null(MR_MED)) rep(NA_real_, K) else as.numeric(MR_MED@Estimate)
    Median.se  <- if (is.null(MR_MED)) rep(NA_real_, K) else as.numeric(MR_MED@StdError)
  }

  EST <- cbind(
    BEE_RB.est, cML_RB.est, MRcML_RB.est,
    BEE.est, MRcML.est, IVW.est, Egger.est, Median.est
  )
  SE <- cbind(
    BEE_RB.se, cML_RB.se, MRcML_RB.se,
    BEE.se, MRcML.se, IVW.se, Egger.se, Median.se
  )
  colnames(EST) <- colnames(SE) <- c("BEE_RB","cML_RB","MRcML_RB","BEE","MRcML","IVW","Egger","Median")

  list(EST = EST, SE = SE, setting = setting, runtime = run.time)
}

## -------------------------
## Parallel plan
## -------------------------
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(ncores) || ncores < 1) ncores <- 1
future::plan(future::multisession, workers = ncores)
cat("Using", ncores, "workers via multisession\n")

t_all <- proc.time()[3]

res_all <- future.apply::future_lapply(
  X = simulation.ind.set,
  FUN = run_one_sim,
  future.seed = TRUE,
  future.packages = c("MASS","MRBEE","MVMRcML","MendelianRandomization","RBCorrection"),
  M = M, K = K, pi1 = pi1, pi2 = pi2, theta = theta,
  Sigma_b_valid = Sigma_b_valid,
  Sigma_b_inval = Sigma_b_inval,
  sigma2y_td = sigma2y_td,
  Sigma_polygenic = Sigma_polygenic,
  Sigma_delta = Sigma_delta,
  nx = nx, ny = ny,
  pthr = pthr, pthr2 = pthr2,
  Thetaxx = Thetaxx,
  Vx = Vx, Vx_inv = Vx_inv,
  eta = eta,
  h2_gamma_perX = h2_gamma_perX,
  h2_xi_target = h2_xi_target,
  RB_HELPER = RB_HELPER
)

total_time <- proc.time()[3] - t_all
cat("Total wall time:", total_time, "seconds\n")

res_list <- lapply(res_all, function(x) list(EST = x$EST, SE = x$SE))
setting  <- lapply(res_all, `[[`, "setting")
runningtime <- lapply(res_all, `[[`, "runtime")

save_file <- paste(
  save_datdir, "simres_MRMethods", set.ind,
  "N", N, "pthr", pthr,
  "theta", paste0(theta, collapse = ","),
  "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio,
  "rhoxy", rhoxy, "h2u", h2x, "h2v", h2y, "rhoos", rhoos, "rhoosx", rhoos_x,
  "h2gamma", h2_gamma_perX, "h2xi", h2_xi_target,
  ".Rdata", sep = ""
)

save(
  res_list,
  setting,
  runningtime,
  total_time,
  simulation.ind.set,
  file = save_file
)

cat("Saved to:", save_file, "\n")