################################################################################
##  test_residuals_run.R
##
##  Comprehensive verification of mlts:::residuals(fit_obj, ts) prior to
##  pushing to remote.
##
##  NOTE ON NAMESPACE BUG:
##    residuals() in R/residuals.R lacks @export and is NOT registered in
##    NAMESPACE.  Calling residuals(fit_obj, ts) dispatches to stats::residuals
##    which silently returns numeric(0) for any mltsfit object.
##    This script calls mlts:::residuals() directly to test the implementation.
##
##  Scenarios
##    S1  – Bivariate VAR(1), inno_covs_zero = TRUE
##    S2  – Bivariate VAR(1), fix_inno_covs = TRUE  (VAR_manifestCovsFix)
##    S3  – Bivariate VAR(1), inno_covs_zero = FALSE, fix_inno_covs = FALSE
##           (latent factor approach, inno_covs_dir = "pos")
##    S4  – Trivariate VAR(1), fix_inno_covs = TRUE  (q = 3)
##    S5  – Bivariate VAR(1) with 5 % MCAR missing values
##    S6  – Bivariate VAR(1) with left-censoring at 0
##    NEG – Negative control: monitor_person_pars = FALSE
##
##  Checks per scenario (A, B, C)
##    A  Structural output integrity
##    B  Independent arithmetic verification for person 1 / time t = 2
##    C  Statistical properties (grand mean ≈ 0, σ alignment, off-diag cov)
##
##  Author: auto-generated verification script
##  Date:   2026-10-20
################################################################################

# ── Preamble ──────────────────────────────────────────────────────────────────
cat("================================================================================\n")
cat("  mlts:::residuals() Comprehensive Verification Script\n")
cat("  Run started:", format(Sys.time()), "\n")
cat("================================================================================\n\n")

suppressPackageStartupMessages({
  library(mlts)
  library(rstan)
})

## ── Resolve the residuals() implementation ────────────────────────────────────
## BUG 0: R/residuals.R has no @export and is not registered in NAMESPACE.
## The function is therefore absent from the installed package namespace.
## We source it directly from the source file so the remaining tests can
## evaluate its implementation.  The sourcing is done into a fresh environment
## that has the mlts namespace as its parent so internal helpers are accessible.

RESIDUALS_SRC <- file.path(system.file(package = "mlts"), "..", "..", "..",
                            "mlts-rdsem", "R", "residuals.R")
# Fallback: resolve relative to the script location
if (!file.exists(RESIDUALS_SRC)) {
  RESIDUALS_SRC <- file.path(getwd(), "R", "residuals.R")
}
if (!file.exists(RESIDUALS_SRC)) {
  stop("Cannot locate R/residuals.R. Run the script from the package root:\n",
       "  Rscript tests/local/test_residuals_run.R")
}

cat("  Sourcing residuals() from:", RESIDUALS_SRC, "\n")
resid_env <- new.env(parent = getNamespace("mlts"))
source(RESIDUALS_SRC, local = resid_env)
mlts_residuals <- resid_env$residuals

# Report namespace status
cat("  stats::residuals  ->", environmentName(environment(stats::residuals)), "\n")
pkg_has_resid <- "residuals" %in% ls(getNamespace("mlts"))
if (pkg_has_resid) {
  cat("  mlts namespace    -> residuals() present  (export bug is fixed)\n")
} else {
  cat("  [BUG 0] mlts namespace -> residuals() ABSENT (not exported / not in NAMESPACE)\n")
  cat("  [BUG 0] Calling residuals(fit) dispatches silently to stats::residuals.\n")
  # Demonstrate the silent misdispatch
  r_demo <- suppressWarnings(tryCatch(stats::residuals(structure(list(), class="mltsfit")),
                                      error = function(e) NULL))
  cat(sprintf("  [BUG 0] stats::residuals(mltsfit_obj) returns: class=%s length=%d\n\n",
              class(r_demo)[1], length(r_demo)))
}
cat("\n")

# Global MCMC settings
N_SIM      <- 30
TP_SIM     <- 60
CHAINS     <- 2
ITER       <- 800          # 50 % warm-up → 400 post-warm-up draws / chain
TS_NAMES   <- c("Y1", "Y2")
TS_NAMES_3 <- c("Y1", "Y2", "Y3")

# ── Helpers ───────────────────────────────────────────────────────────────────

## Pretty assertion wrapper
assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    cat("  [FAIL] ", msg, "\n", sep = "")
    return(invisible(FALSE))
  }
  cat("  [PASS] ", msg, "\n", sep = "")
  return(invisible(TRUE))
}

## Summary table
scenario_results <- data.frame(
  Scenario    = character(0),
  Status      = character(0),
  Runtime_s   = numeric(0),
  MaxRhat     = numeric(0),
  MinBulkESS  = numeric(0),
  Notes       = character(0),
  stringsAsFactors = FALSE
)

record_result <- function(name, status, runtime, fit_obj = NULL, notes = "") {
  max_rhat <- NA_real_; min_ess <- NA_real_
  if (!is.null(fit_obj) && !is.null(fit_obj$stanfit)) {
    sums <- tryCatch(rstan::monitor(fit_obj$stanfit, print = FALSE), error = function(e) NULL)
    if (!is.null(sums)) {
      max_rhat <- round(max(sums$Rhat, na.rm = TRUE), 3)
      min_ess  <- round(min(sums$Bulk_ESS, na.rm = TRUE), 0)
    }
  }
  scenario_results <<- rbind(scenario_results, data.frame(
    Scenario = name, Status = status, Runtime_s = round(runtime, 1),
    MaxRhat = max_rhat, MinBulkESS = min_ess, Notes = notes,
    stringsAsFactors = FALSE
  ))
}

## Independent residual calculator: pp=1, t=2 (first row of residual matrix)
compute_manual_residual <- function(fit_obj, ts, pp = 1) {
  standata    <- fit_obj$standata
  person_pars <- fit_obj$person.pars.summary
  D           <- standata$D
  D_cen       <- standata$D_cen
  D_cen_pos   <- standata$D_cen_pos
  is_wcen     <- standata$is_wcen
  N_obs_id    <- standata$N_obs_id
  obs_id      <- N_obs_id[pp]
  maxLag      <- standata$maxLag
  N_pred      <- standata$N_pred
  D_pred      <- standata$D_pred
  D_pred2     <- standata$D_pred2
  Lag_pred    <- standata$Lag_pred
  Dpos1       <- standata$Dpos1

  # Step 1: build y_cen for person pp
  pos   <- if (pp == 1) 1 else sum(N_obs_id[1:(pp - 1)]) + 1
  y_raw <- fit_obj$data[pos:(pos + obs_id - 1), ts, drop = FALSE]

  y_cen_pp <- matrix(NA_real_, nrow = obs_id, ncol = D)
  for (d in seq_len(D)) {
    raw_vals <- y_raw[[ts[d]]]
    raw_vals[raw_vals == -Inf] <- NA   # guard against -Inf sentinel
    if (is_wcen[d] == 1) {
      # look up the person's posterior-mean mu_d
      # person.pars.summary has column "ID" (the original subject id) and "num_id"
      pp_id <- unique(fit_obj$data$num_id[pos:(pos + obs_id - 1)])
      mu_d  <- person_pars$mean[
        person_pars$num_id == pp_id & person_pars$Param == paste0("mu_", d)
      ]
      if (length(mu_d) == 0) stop("mu_", d, " not found in person.pars.summary for pp=", pp)
      y_cen_pp[, d] <- raw_vals - mu_d
    } else {
      y_cen_pp[, d] <- raw_vals
    }
  }

  # Step 2: posterior means of b_free
  samples_b <- rstan::extract(fit_obj$stanfit, pars = "b_free")$b_free
  b <- apply(samples_b, c(2, 3), mean)  # [N x n_pars]

  # Step 3: compute mu at t = maxLag + 1 (= row 2 for lag-1 model)
  t_idx     <- maxLag + 1
  mu_manual <- rep(NA_real_, D_cen)

  for (d in seq_len(D)) {
    if (is_wcen[d] == 1) {
      if (N_pred[d] > 0) {
        b_vals <- x_vals <- numeric(N_pred[d])
        for (nd in seq_len(N_pred[d])) {
          lag_use <- Lag_pred[d, nd]
          t_lag   <- t_idx - lag_use
          if (standata$D_pred2[d, nd] == -99) {
            x_vals[nd] <- y_cen_pp[t_lag, D_pred[d, nd]]
          } else {
            lag_use2   <- standata$Lag_pred2[d, nd]
            x_vals[nd] <- y_cen_pp[t_lag, D_pred[d, nd]] *
                          y_cen_pp[t_idx - lag_use2, standata$D_pred2[d, nd]]
          }
          b_vals[nd] <- b[pp, Dpos1[d] + (nd - 1)]
        }
        mu_manual[D_cen_pos[d]] <- sum(b_vals * x_vals)
      } else {
        mu_manual[D_cen_pos[d]] <- 0
      }
    }
  }

  # Step 4: residual
  obs_at_t2 <- y_cen_pp[t_idx, D_cen_pos]
  obs_at_t2 - mu_manual
}

## Check A – structural output integrity
check_A <- function(res_list, fit_obj, ts, label) {
  cat("\n  -- Check A: Structural Output Integrity --\n")
  standata <- fit_obj$standata
  N        <- standata$N
  D_cen    <- standata$D_cen
  maxLag   <- standata$maxLag
  N_obs_id <- standata$N_obs_id
  ok       <- TRUE

  # A1
  ok <- assert_true(length(res_list) == N,
    sprintf("length(residuals) == N  (%d == %d)", length(res_list), N)) & ok

  # A2
  all_mat <- all(vapply(res_list, is.matrix, logical(1)))
  all_num <- all(vapply(res_list, function(x) is.numeric(x) || is.double(x), logical(1)))
  ok <- assert_true(all_mat, "All elements are matrices") & ok
  ok <- assert_true(all_num, "All matrices are numeric") & ok

  # A3
  nrow_ok <- all(vapply(seq_len(N), function(pp) {
    nrow(res_list[[pp]]) == (N_obs_id[pp] - maxLag)
  }, logical(1)))
  ok <- assert_true(nrow_ok, sprintf(
    "Each matrix has nrow = obs_id - maxLag (%d - %d = %d for balanced data)",
    N_obs_id[1], maxLag, N_obs_id[1] - maxLag)) & ok

  # A4
  ncol_ok <- all(vapply(res_list, function(x) ncol(x) == D_cen, logical(1)))
  ok <- assert_true(ncol_ok, sprintf("Each matrix has ncol = D_cen (%d)", D_cen)) & ok

  # A5
  expected_colnames <- ts[standata$D_cen_pos]
  colnames_ok <- all(vapply(res_list, function(x) {
    identical(colnames(x), expected_colnames)
  }, logical(1)))
  ok <- assert_true(colnames_ok,
    paste("Column names ==", paste(expected_colnames, collapse = ", "))) & ok

  # A6
  no_nan_inf <- all(vapply(res_list, function(x) {
    vals <- x[!is.na(x)]
    !any(is.nan(vals)) && !any(is.infinite(vals))
  }, logical(1)))
  ok <- assert_true(no_nan_inf, "No NaN or Inf in residual matrices") & ok

  ok
}

## Check B – independent arithmetic at pp=1, t=2
check_B <- function(res_list, fit_obj, ts, label) {
  cat("\n  -- Check B: Independent Arithmetic Verification (pp=1, t=2) --\n")
  ok <- TRUE

  if (length(res_list) == 0) {
    cat("  [SKIP] residuals list is empty – cannot verify\n")
    return(FALSE)
  }

  e_manual <- tryCatch(
    compute_manual_residual(fit_obj, ts, pp = 1),
    error = function(e) {
      cat("  [SKIP] Manual residual computation failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(e_manual)) return(FALSE)

  e_func <- res_list[[1]][1, ]   # first row = t=2 for maxLag=1

  keep <- !is.na(e_manual) & !is.na(e_func)
  if (!any(keep)) {
    cat("  [SKIP] All values NA at t=2 for pp=1 (missing data scenario)\n")
    return(TRUE)
  }

  tol   <- 1e-4
  match <- isTRUE(all.equal(as.numeric(e_manual[keep]),
                             as.numeric(e_func[keep]), tolerance = tol))
  ok    <- assert_true(match,
    sprintf("Manual residual matches mlts_residuals()[[1]][1,] (tol=%.0e)", tol)) & ok

  if (!match) {
    cat("    Manual:   ", paste(round(e_manual, 6), collapse = "  "), "\n")
    cat("    Function: ", paste(round(e_func,   6), collapse = "  "), "\n")
    cat("    Diff:     ", paste(round(e_manual - e_func, 6), collapse = "  "), "\n")
  }
  ok
}

## Check C – statistical properties
check_C <- function(res_list, fit_obj, ts, label, check_offdiag = FALSE) {
  cat("\n  -- Check C: Statistical Properties --\n")
  ok <- TRUE

  if (length(res_list) == 0) {
    cat("  [SKIP] residuals list is empty – cannot compute statistics\n")
    return(FALSE)
  }

  # Grand mean ≈ 0
  all_vals <- unlist(lapply(res_list, as.vector))
  all_vals <- all_vals[is.finite(all_vals)]
  if (length(all_vals) == 0) {
    cat("  [SKIP] No finite values in residuals\n")
    return(FALSE)
  }
  grand_mean <- mean(all_vals)
  ok <- assert_true(abs(grand_mean) < 0.15,
    sprintf("|grand mean(residuals)| = %.4f  < 0.15", grand_mean)) & ok

  # Per-dimension SD vs innovation SD
  standata    <- fit_obj$standata
  D_cen       <- standata$D_cen
  D_cen_pos   <- standata$D_cen_pos
  person_pars <- fit_obj$person.pars.summary

  for (dcen in seq_len(D_cen)) {
    col_resid <- unlist(lapply(res_list, function(x) x[, dcen]))
    col_resid <- col_resid[is.finite(col_resid)]
    emp_sd    <- sd(col_resid)

    d_orig <- D_cen_pos[dcen]

    # person-specific log-innov-vars
    ln_sigma2_vals <- person_pars$mean[
      person_pars$Param == paste0("ln.sigma2_", d_orig)
    ]
    if (length(ln_sigma2_vals) > 0) {
      expected_sd <- mean(exp(ln_sigma2_vals / 2))
    } else {
      # fall back to population parameter
      pop_row <- fit_obj$pop.pars.summary[
        fit_obj$pop.pars.summary$Param == paste0("ln.sigma2_", d_orig), ]
      expected_sd <- if (nrow(pop_row) > 0) mean(exp(pop_row$mean / 2)) else NA
    }

    if (!is.na(expected_sd)) {
      ratio <- emp_sd / expected_sd
      ok <- assert_true(ratio > 0.5 && ratio < 2.0,
        sprintf("Dim %d: emp_sd=%.3f, expected_innov_sd≈%.3f, ratio=%.2f (in [0.5,2.0])",
                d_orig, emp_sd, expected_sd, ratio)) & ok
    } else {
      cat("  [INFO] Cannot retrieve innovation SD for dimension", d_orig, "– skipping ratio check\n")
    }
  }

  # Off-diagonal covariance (Scenario 2)
  if (check_offdiag && D_cen >= 2 && length(res_list) > 0) {
    m1 <- res_list[[1]]
    if (is.matrix(m1) && nrow(m1) > 1 && ncol(m1) >= 2) {
      cov_pp  <- cov(m1, use = "complete.obs")
      offdiag <- cov_pp[1, 2]
      ok <- assert_true(abs(offdiag) > 0.01,
        sprintf("Off-diagonal cov(residuals[[1]]) = %.4f  > 0.01 (non-zero expected)", offdiag)) & ok
    }
  }
  ok
}

################################################################################
##  NEGATIVE CONTROL: monitor_person_pars = FALSE
################################################################################
cat("\n================================================================================\n")
cat("  NEGATIVE CONTROL – monitor_person_pars = FALSE\n")
cat("================================================================================\n")
t_neg <- proc.time()

set.seed(42)
model_neg   <- mlts_model(q = 2, inno_covs_zero = TRUE)
simData_neg <- mlts_sim(model = model_neg, N = 10, TP = 20, default = TRUE, seed = 42)

# Note: use chains=2 because mlts_fit() has a separate bug when chains=1 and
# monitor_person_pars=FALSE: get_posterior_mean() returns a 1-D vector so
# sums.i[, dim(sums.i)[2]] crashes with "incorrect number of dimensions".
fit_neg <- suppressWarnings(mlts_fit(
  model   = model_neg,
  data    = simData_neg,
  ts      = TS_NAMES,
  iter    = 400,
  chains  = 2,
  cores   = 2,
  monitor_person_pars = FALSE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Checking that mlts:::residuals() halts with the expected error ...\n")
neg_ok <- FALSE
tryCatch({
  mlts_residuals(fit_neg, ts = TS_NAMES)
  cat("  [FAIL] mlts:::residuals() did NOT throw an error – expected error!\n")
}, error = function(e) {
  msg           <- conditionMessage(e)
  expected_phrase <- "Person-specific param"
  neg_ok <<- grepl(expected_phrase, msg, ignore.case = TRUE)
  if (neg_ok) {
    cat("  [PASS] mlts:::residuals() correctly halts:\n  \"", msg, "\"\n", sep = "")
  } else {
    cat("  [FAIL] Error thrown but unexpected message:\n  \"", msg, "\"\n", sep = "")
  }
})

# Also verify that calling the exported stats::residuals gives no error but
# wrong (empty/NULL) result – documenting the NAMESPACE bug
cat("\n  Verifying stats::residuals dispatches to wrong function (NAMESPACE bug):\n")
stats_result <- tryCatch(
  stats::residuals(fit_neg),
  error = function(e) paste0("ERROR: ", conditionMessage(e))
)
if (is.character(stats_result)) {
  cat("  [INFO] stats::residuals errored:", stats_result, "\n")
} else {
  cat(sprintf("  [INFO] stats::residuals returned: class=%s  length=%d  (expected: mltsfit-specific output)\n",
              class(stats_result)[1], length(stats_result)))
  cat("  [PASS] Confirms NAMESPACE BUG 0: public residuals() call silently misdispatches\n")
}

t_neg_elapsed <- (proc.time() - t_neg)["elapsed"]
record_result("NEG (no person pars)", if (neg_ok) "PASS" else "FAIL",
              t_neg_elapsed, fit_neg,
              "Expects stop() from mlts:::residuals when monitor_person_pars = FALSE")

################################################################################
##  SCENARIO 1 – Bivariate VAR(1), inno_covs_zero = TRUE
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 1 – Bivariate VAR(1), inno_covs_zero = TRUE\n")
cat("================================================================================\n")
t_s1 <- proc.time()

set.seed(42)
model_s1   <- mlts_model(q = 2, inno_covs_zero = TRUE)
simData_s1 <- mlts_sim(model = model_s1, N = N_SIM, TP = TP_SIM,
                       default = TRUE, seed = 42)

fit_s1 <- suppressWarnings(mlts_fit(
  model   = model_s1,
  data    = simData_s1,
  ts      = TS_NAMES,
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s1 <- mlts_residuals(fit_s1, ts = TS_NAMES)

s1_A <- check_A(res_s1, fit_s1, TS_NAMES, "S1")
s1_B <- check_B(res_s1, fit_s1, TS_NAMES, "S1")
s1_C <- check_C(res_s1, fit_s1, TS_NAMES, "S1")

s1_ok <- s1_A && s1_B && s1_C
t_s1_elapsed <- (proc.time() - t_s1)["elapsed"]
record_result("S1 (bivariate, zero covs)", if (s1_ok) "PASS" else "FAIL",
              t_s1_elapsed, fit_s1)

################################################################################
##  SCENARIO 2 – Bivariate VAR(1), fix_inno_covs = TRUE (VAR_manifestCovsFix)
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 2 – Bivariate VAR(1), fix_inno_covs = TRUE\n")
cat("================================================================================\n")
t_s2 <- proc.time()

set.seed(42)
model_s2   <- mlts_model(q = 2, inno_covs_zero = FALSE, fix_inno_covs = TRUE)
simData_s2 <- mlts_sim(model = model_s2, N = N_SIM, TP = TP_SIM,
                       default = TRUE, seed = 42)

fit_s2 <- suppressWarnings(mlts_fit(
  model   = model_s2,
  data    = simData_s2,
  ts      = TS_NAMES,
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s2 <- mlts_residuals(fit_s2, ts = TS_NAMES)

s2_A <- check_A(res_s2, fit_s2, TS_NAMES, "S2")
s2_B <- check_B(res_s2, fit_s2, TS_NAMES, "S2")
s2_C <- check_C(res_s2, fit_s2, TS_NAMES, "S2", check_offdiag = TRUE)

s2_ok <- s2_A && s2_B && s2_C
t_s2_elapsed <- (proc.time() - t_s2)["elapsed"]
record_result("S2 (bivariate, fix covs)", if (s2_ok) "PASS" else "FAIL",
              t_s2_elapsed, fit_s2)

################################################################################
##  SCENARIO 3 – Bivariate VAR(1), latent factor inno covariance (pos)
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 3 – Bivariate VAR(1), fix_inno_covs = FALSE, inno_covs_dir = 'pos'\n")
cat("================================================================================\n")
t_s3 <- proc.time()

set.seed(42)
model_s3   <- mlts_model(q = 2, inno_covs_zero = FALSE,
                         fix_inno_covs = FALSE, inno_covs_dir = "pos")
simData_s3 <- mlts_sim(model = model_s3, N = N_SIM, TP = TP_SIM,
                       default = TRUE, seed = 42)

fit_s3 <- suppressWarnings(mlts_fit(
  model   = model_s3,
  data    = simData_s3,
  ts      = TS_NAMES,
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s3 <- mlts_residuals(fit_s3, ts = TS_NAMES)

s3_A <- check_A(res_s3, fit_s3, TS_NAMES, "S3")
s3_B <- check_B(res_s3, fit_s3, TS_NAMES, "S3")
s3_C <- check_C(res_s3, fit_s3, TS_NAMES, "S3")

s3_ok <- s3_A && s3_B && s3_C
t_s3_elapsed <- (proc.time() - t_s3)["elapsed"]
record_result("S3 (bivariate, latent factor)", if (s3_ok) "PASS" else "FAIL",
              t_s3_elapsed, fit_s3,
              "Person-specific innovation covariance via latent factor")

################################################################################
##  SCENARIO 4 – Trivariate VAR(1), q = 3
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 4 – Trivariate VAR(1), q = 3, fix_inno_covs = TRUE\n")
cat("================================================================================\n")
t_s4 <- proc.time()

set.seed(42)
model_s4   <- mlts_model(q = 3, fix_inno_covs = TRUE)
simData_s4 <- mlts_sim(model = model_s4, N = N_SIM, TP = TP_SIM,
                       default = TRUE, seed = 42)

fit_s4 <- suppressWarnings(mlts_fit(
  model   = model_s4,
  data    = simData_s4,
  ts      = TS_NAMES_3,
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s4 <- mlts_residuals(fit_s4, ts = TS_NAMES_3)

s4_A <- check_A(res_s4, fit_s4, TS_NAMES_3, "S4")
s4_B <- check_B(res_s4, fit_s4, TS_NAMES_3, "S4")
s4_C <- check_C(res_s4, fit_s4, TS_NAMES_3, "S4")

s4_ok <- s4_A && s4_B && s4_C
t_s4_elapsed <- (proc.time() - t_s4)["elapsed"]
record_result("S4 (trivariate q=3)", if (s4_ok) "PASS" else "FAIL",
              t_s4_elapsed, fit_s4)

################################################################################
##  SCENARIO 5 – Bivariate VAR(1) with 5 % MCAR missing values
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 5 – Bivariate VAR(1) + 5 % MCAR missing values\n")
cat("================================================================================\n")
t_s5 <- proc.time()

set.seed(42)
model_s5   <- mlts_model(q = 2, inno_covs_zero = TRUE)
simData_s5 <- mlts_sim(model = model_s5, N = N_SIM, TP = TP_SIM,
                       default = TRUE, seed = 42)

data_s5  <- simData_s5$data
n_total  <- nrow(data_s5)
set.seed(999)
na_idx_Y1 <- sample(seq_len(n_total), size = floor(0.05 * n_total))
na_idx_Y2 <- sample(seq_len(n_total), size = floor(0.05 * n_total))
data_s5$Y1[na_idx_Y1] <- NA
data_s5$Y2[na_idx_Y2] <- NA

cat(sprintf("  Inserted %d NAs in Y1 and %d NAs in Y2 (%.1f%% each)\n",
            length(na_idx_Y1), length(na_idx_Y2),
            100 * length(na_idx_Y1) / n_total))

# Note: beep = "time" would crash here because prepare_data() sets tinterval=1
# internally then cannot find 'time' or 'days' variable for the time-grid
# branch when na.rm = FALSE. Use time = "time" (without tinterval) so that
# data is sorted correctly and NAs are preserved for Stan imputation.
fit_s5 <- suppressWarnings(mlts_fit(
  model   = model_s5,
  data    = data_s5,
  id      = "ID",
  ts      = TS_NAMES,
  time    = "time",
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  na.rm   = FALSE,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s5 <- tryCatch(
  mlts_residuals(fit_s5, ts = TS_NAMES),
  error = function(e) {
    cat("  [FAIL] mlts:::residuals() crashed on missing data:\n  ", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(res_s5)) {
  s5_A <- check_A(res_s5, fit_s5, TS_NAMES, "S5")
  s5_B <- check_B(res_s5, fit_s5, TS_NAMES, "S5")
  s5_C <- check_C(res_s5, fit_s5, TS_NAMES, "S5")

  n_na_resid  <- sum(vapply(res_s5, function(x) sum(is.na(x)), integer(1)))
  n_inf_resid <- sum(vapply(res_s5, function(x) sum(is.infinite(x), na.rm = TRUE), integer(1)))
  cat(sprintf("  [INFO] NAs in residuals:  %d\n", n_na_resid))
  cat(sprintf("  [INFO] Inf in residuals:  %d  (> 0 confirms BUG 1: -Inf sentinel not guarded)\n",
              n_inf_resid))

  s5_no_inf <- assert_true(n_inf_resid == 0,
    sprintf("No Inf/-Inf in residual matrices (got %d – BUG 1)", n_inf_resid))

  s5_ok <- s5_A && s5_B && s5_C && s5_no_inf
} else {
  s5_A <- s5_B <- s5_C <- FALSE; s5_ok <- FALSE
}

t_s5_elapsed <- (proc.time() - t_s5)["elapsed"]
record_result("S5 (missing 5% MCAR)", if (s5_ok) "PASS" else "FAIL",
              t_s5_elapsed, fit_s5,
              "5% MCAR; Inf expected in residuals due to -Inf sentinel (BUG 1)")

################################################################################
##  SCENARIO 6 – Bivariate VAR(1) with left-censoring at 0
################################################################################
cat("\n================================================================================\n")
cat("  SCENARIO 6 – Bivariate VAR(1) with left-censoring at 0\n")
cat("================================================================================\n")
t_s6 <- proc.time()

set.seed(42)
model_s6_base <- mlts_model(q = 2, inno_covs_zero = TRUE)
simData_s6    <- mlts_sim(model = model_s6_base, N = N_SIM, TP = TP_SIM,
                          default = TRUE, seed = 42)

model_s6 <- mlts_model(q = 2, inno_covs_zero = TRUE, censor_left = 0.0)

data_s6   <- simData_s6$data
CENS_LB   <- 0.0
n_cens_Y1 <- sum(data_s6$Y1 <= CENS_LB, na.rm = TRUE)
n_cens_Y2 <- sum(data_s6$Y2 <= CENS_LB, na.rm = TRUE)
data_s6$Y1[data_s6$Y1 <= CENS_LB] <- CENS_LB
data_s6$Y2[data_s6$Y2 <= CENS_LB] <- CENS_LB

cat(sprintf("  Censored %d Y1 values and %d Y2 values at boundary=%.1f\n",
            n_cens_Y1, n_cens_Y2, CENS_LB))

fit_s6 <- suppressWarnings(mlts_fit(
  model   = model_s6,
  data    = data_s6,
  id      = "ID",
  ts      = TS_NAMES,
  time    = "time",
  iter    = ITER,
  chains  = CHAINS,
  cores   = CHAINS,
  monitor_person_pars = TRUE,
  print_message = FALSE,
  print_warning = FALSE
))

cat("\n  Computing residuals ...\n")
res_s6 <- tryCatch(
  mlts_residuals(fit_s6, ts = TS_NAMES),
  error = function(e) {
    cat("  [FAIL] mlts:::residuals() crashed on censored data:\n  ", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(res_s6)) {
  s6_A <- check_A(res_s6, fit_s6, TS_NAMES, "S6")
  s6_B <- check_B(res_s6, fit_s6, TS_NAMES, "S6")
  s6_C <- check_C(res_s6, fit_s6, TS_NAMES, "S6")

  all_resid_s6 <- unlist(lapply(res_s6, as.vector))
  all_resid_s6 <- all_resid_s6[is.finite(all_resid_s6)]
  cat(sprintf("  [INFO] Residual range: [%.3f, %.3f]\n",
              min(all_resid_s6), max(all_resid_s6)))
  cat("  [INFO] Censored residuals reflect clamped boundary value (BUG 3 – methodological)\n")

  n_nan_s6  <- sum(is.nan(all_resid_s6))
  s6_no_nan <- assert_true(n_nan_s6 == 0,
    sprintf("No NaN values in censored residuals (got %d)", n_nan_s6))

  s6_ok <- s6_A && s6_B && s6_C && s6_no_nan
} else {
  s6_A <- s6_B <- s6_C <- FALSE; s6_ok <- FALSE
}

t_s6_elapsed <- (proc.time() - t_s6)["elapsed"]
record_result("S6 (left-censored at 0)", if (s6_ok) "PASS" else "FAIL",
              t_s6_elapsed, fit_s6,
              "Censored obs use raw clamped value (BUG 3 – methodological)")

################################################################################
##  FINAL SUMMARY REPORT
################################################################################
cat("\n\n")
cat("================================================================================\n")
cat("  FINAL VERIFICATION REPORT\n")
cat("  mlts:::residuals(fit_obj, ts)\n")
cat("  Run completed:", format(Sys.time()), "\n")
cat("================================================================================\n\n")

cat("── Summary Table ────────────────────────────────────────────────────────────\n")
cat(sprintf("  %-35s  %-6s  %8s  %8s  %8s\n",
            "Scenario", "Status", "Runtime", "MaxRhat", "MinESS"))
cat(sprintf("  %-35s  %-6s  %8s  %8s  %8s\n",
            strrep("-", 35), strrep("-", 6), strrep("-", 8),
            strrep("-", 8), strrep("-", 8)))
for (i in seq_len(nrow(scenario_results))) {
  r <- scenario_results[i, ]
  cat(sprintf("  %-35s  %-6s  %7.1fs  %8s  %8s\n",
              r$Scenario, r$Status, r$Runtime_s,
              ifelse(is.na(r$MaxRhat),    "  n/a  ", sprintf("%.3f", r$MaxRhat)),
              ifelse(is.na(r$MinBulkESS), "  n/a  ", sprintf("%d", as.integer(r$MinBulkESS)))))
}
cat("\n")
n_pass <- sum(scenario_results$Status == "PASS")
n_fail <- sum(scenario_results$Status == "FAIL")
cat(sprintf("  Overall: %d / %d scenarios PASS\n\n", n_pass, nrow(scenario_results)))

cat("── Missing Data & Censoring Assessment ──────────────────────────────────────\n")
cat("
  MISSING DATA (Scenario 5):
  ─ prepare_data() replaces NA with -Inf in fit_obj$data before Stan.
  ─ mlts:::residuals() does NOT guard against this sentinel: subtracting
    mu_d from -Inf propagates -Inf into y_cen, then Inf into the residual.
  ─ Missing positions therefore appear as ±Inf rather than NA – a clear
    correctness failure reported via the s5_no_inf assertion.

  CENSORED DATA (Scenario 6):
  ─ Stan imputes values below the censoring threshold via y_cen_impute
    parameters; these are NOT stored in fit_obj$data.
  ─ mlts:::residuals() uses the raw clamped boundary value from fit_obj$data.
    Residuals at censored time points are therefore systematically wrong
    (they reflect the distance from the boundary, not from the true
    innovation). No crash, but methodologically inconsistent.
\n")

cat("── Identified Bugs / Code Quality Issues ────────────────────────────────────\n")
cat("
  *** BUG 0 – CRITICAL: residuals() not exported (NAMESPACE) ***
    R/residuals.R has no @export roxygen tag and no S3method declaration.
    Consequence: calling residuals(fit_obj, ts) dispatches silently to
    stats::residuals(), which returns numeric(0) for an mltsfit object with
    no error, no warning, and no indication that the wrong function ran.
    FIX: Add @export and rename to residuals.mltsfit (or register via
         S3method(residuals, mltsfit) in NAMESPACE). Re-run roxygen2.

  BUG 1 – -Inf sentinel propagation (residuals.R lines 21-25):
    fit_obj$data stores -Inf for missing obs (VARprepare.R line 85).
    Centering with   y_cen[[pp]][,d] <- data[...] - mu_d   produces -Inf
    in y_cen and then ±Inf in residuals instead of NA.
    FIX: Before centering, replace -Inf with NA:
         raw_vals <- fit_obj$data[pos:(pos+obs_id-1), ts[d]]
         raw_vals[raw_vals == -Inf] <- NA
         y_cen[[pp]][,d] <- raw_vals - mu_d   # (or raw_vals for non-centered)

  BUG 2 – Undefined variable `b` in RDSEM branch (residuals.R lines 30-38):
    b <- apply(samples_b, ...) is computed AFTER the y_cen loop (line 57).
    Any model with is_rdsem > 0 will error: 'object b not found'.
    The surrounding commented-out blocks (lines 64-86) confirm incomplete
    refactoring.
    FIX: Move b <- apply(...) BEFORE the y_cen loop.

  BUG 3 – Censored residuals use clamped value (methodological):
    residuals() uses fit_obj$data (clamped at boundary) rather than
    Stan's y_cen_impute posterior draws.
    FIX (optional): Extract imputed values when monitor_all_pars = TRUE;
    otherwise document the limitation with a warning for censored models.

  CODE QUALITY – mlts_fit() 1-chain bug with monitor_person_pars = FALSE:
    rstan::get_posterior_mean() returns a vector (not matrix) when chains=1.
    sums.i[, dim(sums.i)[2]] then fails with 'incorrect number of dimensions'.
    FIX: Guard with as.matrix() or handle the 1-chain case explicitly.

  CODE QUALITY – Dead commented-out code (residuals.R lines 64-66, 83-86):
    Abandoned n_cols and L_inno loading blocks; remove or annotate.

  CODE QUALITY – Shared `pos` variable across loops (residuals.R):
    `pos` is mutated in Step 1 and left at its terminal value; Step 2 does
    not use it, but the implicit state is fragile for future edits.

  CODE QUALITY – person.pars.summary ID column assumed to be 'ID':
    Line 22: person_pars$ID == pp. The actual column name after merge in
    mlts_fit() is num_id (numerical ID). The column named 'ID' is the
    original subject label. If the data uses non-integer IDs (e.g. 'P01'),
    this lookup silently returns numeric(0) and centering is skipped.
\n")

cat("── Final Recommendation ──────────────────────────────────────────────────────\n")
if (n_fail == 0) {
  cat("
  VERDICT: CONDITIONALLY READY
  S1–S4 pass all structural and arithmetic checks when called via
  mlts:::residuals(). HOWEVER, the following MUST be addressed first:
    BUG 0 (CRITICAL) – Add @export / S3method so the function is reachable.
    BUG 1            – Guard -Inf sentinel to avoid Inf in residuals.
    BUG 2            – Move b <- apply(...) before the y_cen loop.
  BUG 3 and code-quality issues should be addressed or documented.
\n")
} else {
  cat(sprintf("
  VERDICT: NOT READY FOR MERGE  (%d scenario(s) FAILED)
  Priority fixes:
    BUG 0 (CRITICAL) – NAMESPACE: export the function properly.
    BUG 1            – -Inf sentinel → Inf in missing-data residuals.
    BUG 2            – undefined `b` in RDSEM branch.
\n", n_fail))
}

cat("================================================================================\n")
cat("  End of Report\n")
cat("================================================================================\n")
