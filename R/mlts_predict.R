
#' @title Predictions for Multilevel Latent Time Series Models
#'
#' @description
#' Calculate fitted conditional means on the original observation scale from a
#' fitted \code{mlts} model. Return posterior summaries or predictions calculated
#' separately for each posterior draw.
#'
#' @param fit_obj A fitted model object of class \code{mlts.fit}, as returned by
#'   \code{\link{mlts_fit}}. Required; no default.
#' @param ts Character vector of time-series variable names in
#'   \code{fit_obj$data}, in the order used by the fitted model, including
#'   time-varying predictors. Required; no default. Output labels use the first
#'   \code{D_cen} entries, where \code{D_cen = fit_obj$standata$D_cen}.
#' @param summary Logical. If \code{TRUE} (default), return posterior mean
#'   predictions and optionally posterior standard errors and credible intervals.
#'   If \code{FALSE}, return arrays containing predictions for individual draws.
#' @param se Logical. If \code{TRUE} (default), include posterior standard
#'   deviations in the summarized \code{"data.frame"} output. Ignored when
#'   \code{summary = FALSE}; disabled when person-specific draws are unavailable.
#' @param ci Logical. If \code{TRUE} (default), include equal-tailed credible
#'   intervals in the summarized \code{"data.frame"} output. Ignored when
#'   \code{summary = FALSE}; disabled when person-specific draws are unavailable.
#' @param ci_level Numeric scalar between 0 and 1 specifying the credible
#'   interval probability. Defaults to \code{0.95}. Used only when credible
#'   intervals are calculated.
#' @param output Character string specifying the output format:
#'   \code{"data.frame"} (default) or \code{"matrix_list"}. The latter requires
#'   \code{summary = TRUE} and returns point estimates only. With
#'   \code{summary = FALSE}, use \code{"data.frame"}; the returned object is
#'   nevertheless a list of arrays and an index data frame.
#'
#' @details
#' Predictions are calculated for the observations stored in \code{fit_obj$data}.
#' For each posterior draw, modeled outcome variables are centered using their
#' person-specific means, whereas uncentered predictors retain their observed
#' values. The within-person conditional mean is calculated from the model's
#' predictor terms, including their specified lags and product terms. The
#' person-specific mean is then added to obtain predictions on the original
#' observation scale. Fixed parameters are assigned according to group membership.
#'
#' For residual dynamic structural equation models, contemporaneous effects are
#' removed before calculating the dynamic component and added back when
#' constructing predictions on the observation scale. The returned \code{mu}
#' component contains the dynamic conditional mean before adding person-specific
#' means and these contemporaneous effects. An outcome without predictor terms
#' has a dynamic conditional mean of zero.
#'
#' No innovations are simulated. With posterior draws available, the reported
#' SE is the standard deviation of the conditional mean predictions across
#' draws, not a Monte Carlo standard error. Credible limits are the posterior
#' quantiles at \code{(1 - ci_level) / 2} and
#' \code{1 - (1 - ci_level) / 2}. These quantities describe parameter uncertainty
#' in fitted conditional means, rather than posterior predictive uncertainty
#' for new observations.
#'
#' If the model contains random parameters but person-specific \code{b_free}
#' draws are unavailable, as when fitted with
#' \code{monitor_person_pars = FALSE}, a warning is issued. Predictions are then
#' calculated once using posterior mean person-specific parameters and, when
#' present, posterior mean fixed parameters obtained from the fitted Stan object.
#' This requires those posterior means to be available. The resulting prediction
#' at posterior mean parameters need not equal the posterior mean prediction.
#' Because there is no distribution of draw-specific predictions in this mode,
#' posterior SEs and credible intervals cannot be computed and are disabled.
#' With \code{summary = FALSE}, the draw dimension has length one.
#'
#' The first \code{fit_obj$standata$maxLag} observations of every individual are
#' retained as \code{NA}, including for outcomes with shorter individual lag
#' requirements. Rows follow the stored data order; the implementation assumes
#' contiguous individual records consistent with the fitted observation counts.
#' Each individual must have more observations than the maximum lag.
#'
#' Only the \code{D_cen} modeled, within-person-centered outcomes are returned.
#' Latent multi-indicator models are not supported and cause an error. The
#' function does not accept new data or produce forecasts beyond the fitted
#' observations, impute missing observations, or apply censoring adjustments.
#' Missing input values can therefore propagate into predictions.
#' Draw-based uncertainty summaries require multiple posterior draws.
#' When \code{summary = FALSE}, \code{se}, \code{ci}, and \code{ci_level} are
#' ignored; a warning is issued if \code{se} or \code{ci} is \code{TRUE}.
#'
#' @return Let \code{M} be \code{nrow(fit_obj$data)}, \code{N} the number of
#'   individuals, \code{D_cen} the number of modeled outcomes, and \code{S} the
#'   number of extracted posterior draws across chains (one in fallback mode).
#'   The return value depends on \code{summary} and \code{output}:
#' \describe{
#'   \item{\code{summary = TRUE, output = "data.frame"}}{
#'     A data frame with \code{M} rows, columns \code{num_id} and
#'     \code{obs_index}, and one \code{<variable>_pred} column per outcome.
#'     Optional \code{<variable>_se}, \code{<variable>_lower}, and
#'     \code{<variable>_upper} columns contain posterior standard deviations and
#'     credible limits. \code{num_id} is the stored numeric individual ID;
#'     \code{obs_index} counts observations within each individual from one.}
#'   \item{\code{summary = TRUE, output = "matrix_list"}}{
#'     An unnamed list of length \code{N}, in fitted individual order. Element
#'     \code{i} is a numeric matrix with \code{N_obs_id[i]} rows and
#'     \code{D_cen} columns containing point predictions. Matrices have no
#'     variable column names. No SEs, credible limits, or index columns are
#'     returned.}
#'   \item{\code{summary = FALSE}}{
#'     A list with components \code{mu}, \code{y_pred}, and \code{index}.
#'     Both arrays have dimensions \code{S} by \code{M} by \code{D_cen},
#'     corresponding to draws, observations, and outcomes. \code{mu} contains
#'     dynamic conditional means and \code{y_pred} contains predictions on the
#'     original observation scale. Their third dimensions are named using
#'     \code{ts[seq_len(D_cen)]}. \code{index} is an \code{M}-row data frame
#'     with columns \code{ID} and \code{obs_index}, identifying the observations
#'     along the second array dimension.}
#' }
#'
#' @seealso \code{\link{mlts_fit}}, \code{\link{mlts_residuals}},
#'   \code{\link{mlts_posterior_sample}}.
#' @export

mlts_predict <- function(fit_obj, ts, summary = TRUE, se = TRUE, ci = TRUE, ci_level = 0.95, output = "data.frame"){
  std <- fit_obj$standata

  if ("Measurement" %in% fit_obj$model$Model) {
    stop("This function is currently not supported for latent multi-indicator models.", call. = FALSE)
  }

  if (!summary && (se || ci)) {
    warning("Arguments 'ci', 'se' and 'ci_level' will be ignored when summary = FALSE.", call. = FALSE)
  }
  has_b_free <- any(grepl("^b_free", names(fit_obj$stanfit)))
  fallback <- FALSE
  if (std$n_random > 0 && !has_b_free){
    fallback <- TRUE
    se <- FALSE
    ci <- FALSE
    warning("Person-specific posterior draws are not available.
     Predictions and residuals will be calculated using posterior means.
      Posterior standard errors and credible intervals cannot be computed.
       Consider setting 'monitor_person_pars' = TRUE during fitting.")
  }
  if (output != "data.frame" && output != "matrix_list"){
    stop("'output' has to be set to 'data.frame' or 'matrix_list'.", call. = FALSE)
  }
  if (output == "matrix_list" && summary == FALSE) {
    stop("output 'matrix_list' is only possible when summary = TRUE ")
  }


  n_fixed <- std$n_fixed
  is_fixed <- std$is_fixed
  n_random <- std$n_random
  is_random <- std$is_random
  N_G <- std$N_G
  N <- std$N
  P <- std$n_pars - std$n_innos_fix
  G <- std$G
  maxLag <- std$maxLag
  D_cen <- std$D_cen
  N_obs_id <- std$N_obs_id
  D <- std$D

  if (fallback){
    S <- 1
    b_free_col <- utils::tail(colnames(rstan::get_posterior_mean(fit_obj$stanfit, "b_free")), 1)
    b_free_samples <- rstan::get_posterior_mean(fit_obj$stanfit, "b_free")[, b_free_col]

    b_free <- array(NA_real_, dim = c(1, N, n_random))
    b_free[1, , ] <- matrix(b_free_samples, nrow = N, ncol = n_random, byrow = TRUE)

    if (n_fixed > 0){
      b_fix_col <- tail(colnames(rstan::get_posterior_mean(fit_obj$stanfit, "b_fix")), 1)
      b_fix_samples <- rstan::get_posterior_mean(fit_obj$stanfit, "b_fix")[, b_fix_col]
      
      b_fix <- array(NA_real_, dim = c(1, G, n_fixed))
      b_fix[1, , ] <- matrix(b_fix_samples, nrow = G, ncol = n_fixed, byrow = TRUE)
    }
  } else {
    if (n_random > 0) {
      samples_b_free <- rstan::extract(fit_obj$stanfit, pars = "b_free")$b_free
      if (length(dim(samples_b_free)) == 2) {
        dim(samples_b_free) <- c(dim(samples_b_free)[1], N, n_random)
      }
      S_total <- dim(samples_b_free)[1]
      b_free <- samples_b_free
    }

    if (n_fixed > 0){
      samples_b_fix <- rstan::extract(fit_obj$stanfit, pars = "b_fix")$b_fix
      if (length(dim(samples_b_fix)) == 2) {
        dim(samples_b_fix) <- c(dim(samples_b_fix)[1], G, n_fixed)
      }
      S_total <- dim(samples_b_fix)[1]
      b_fix <- samples_b_fix
    }
    S <- S_total
  }

  b <- array(data = NA_real_, dim = c(S, N, P))
  if (n_random > 0){
    b[, , is_random] <- b_free
  }

  if (n_fixed > 0){
    for (s in seq_len(S)){
      for (i in seq_len(n_fixed)){
        for (g in seq_len(G)){
          b[s, std$g_id_pos[g, 1:N_G[g]], is_fixed[1, i]] = b_fix[s, g, i]
        }
      }
    }
  }
  mus_all <- vector("list", N)
  y_pred_all <- vector("list", N)

  for (pp in seq_len(N)){
    obs_p <- N_obs_id[pp]
    mus_all[[pp]] <- array(NA_real_, dim = c(S, obs_p - maxLag, D_cen))
    y_pred_all[[pp]] <- array(NA_real_, dim = c(S, obs_p - maxLag, D_cen))
  }
  for (s in seq_len(S)){
    result <- predict_iteration(fit_obj, ts, b[s, , ])

    for (pp in seq_len(N)){
      mus_all[[pp]][s, , ] <- result$mus[[pp]]
      y_pred_all[[pp]][s, , ] <- result$y_pred[[pp]]
    }
  }
  obs_index <- stats::ave(
      seq_len(nrow(fit_obj$data)),
      fit_obj$data$num_id,
      FUN = seq_along
    )

  index_df <- data.frame(ID = fit_obj$data$num_id,
                         obs_index = obs_index)
  if (summary == TRUE){
    y_pred_means <- vector("list", N)

    if (se) {
      y_pred_ses <- vector("list", N)
      y_pred_mat_se <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    }
    if (ci) {
      y_pred_lower <- vector("list", N)
      y_pred_upper <- vector("list", N)
      pred_mat_lower <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
      pred_mat_upper <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    }
    if (S > 1){
      probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
      for (pp in seq_len(N)){
        y_pred_means[[pp]] <- apply(y_pred_all[[pp]], c(2, 3), mean)
        if (se == TRUE){
          y_pred_ses[[pp]] <- apply(y_pred_all[[pp]], c(2, 3), stats::sd)
        }
        if (ci == TRUE) {
          q_pp <- apply(y_pred_all[[pp]], c(2, 3), stats::quantile, probs = probs, na.rm = TRUE)
          y_pred_lower[[pp]] <- array(q_pp[1, , ], dim = dim(y_pred_all[[pp]])[2:3])
          y_pred_upper[[pp]] <- array(q_pp[2, , ], dim = dim(y_pred_all[[pp]])[2:3])
        }
      }
    } else {
      for (pp in seq_len(N)) {
        y_pred_means[[pp]] <- array(y_pred_all[[pp]][1, , , drop = FALSE], dim = dim(y_pred_all[[pp]])[2:3])
      }
    }

    y_pred_mat <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    pos = 1
    for (pp in seq_len(N)){
      n_obs <- N_obs_id[pp]
      y_pred_mat[(pos + maxLag):(pos + n_obs - 1), ] <- y_pred_means[[pp]]
      if (se) y_pred_mat_se[(pos + maxLag):(pos + n_obs - 1), ] <- y_pred_ses[[pp]]
      if (ci) {
        pred_mat_lower[(pos + maxLag):(pos + n_obs - 1), ] <- y_pred_lower[[pp]]
        pred_mat_upper[(pos + maxLag):(pos + n_obs - 1), ] <- y_pred_upper[[pp]]
      }
      pos <- pos + n_obs
    }

    if (output == "matrix_list"){
      pos <- 1
      result_list <- vector("list", N)
      for (pp in seq_len(N)){
        n_obs <- N_obs_id[pp]
        result_list[[pp]] <- y_pred_mat[pos : (pos + n_obs - 1), , drop = FALSE]
        pos <- pos + n_obs
      }
      return(result_list)
    } else if (output == "data.frame") {
      result_df <- index_df
    colnames(result_df) <- c("num_id", "obs_index")
    y_pred_df <- as.data.frame(y_pred_mat)
    colnames(y_pred_df) <- paste0(ts[seq_len(D_cen)], "_pred")
    result_df <- cbind(result_df, y_pred_df)

    if (se) {
      se_df <- as.data.frame(y_pred_mat_se)
      colnames(se_df) <- paste0(ts[seq_len(D_cen)], "_se")
      result_df <- cbind(result_df, se_df)
    }
    if (ci) {
      lower_df <- as.data.frame(pred_mat_lower)
      colnames(lower_df) <- paste0(ts[seq_len(D_cen)], "_lower")
      upper_df <- as.data.frame(pred_mat_upper)
      colnames(upper_df) <- paste0(ts[seq_len(D_cen)], "_upper")
      result_df <- cbind(result_df, lower_df, upper_df)
    }
    return(result_df)
    }

  } else { 
    out_mat_mu <- array(NA_real_, dim = c(S, nrow(fit_obj$data), D_cen))
    out_mat_y_pred <- array(NA_real_, dim = c(S, nrow(fit_obj$data), D_cen))
    pos = 1
    for (pp in seq_len(N)){
      n_obs <- N_obs_id[pp]
      out_mat_mu[, (pos + maxLag):(pos + n_obs - 1), ] <- mus_all[[pp]]
      out_mat_y_pred[, (pos + maxLag):(pos + n_obs - 1), ] <- y_pred_all[[pp]]
      pos = pos + n_obs
    }
    dimnames(out_mat_mu)[[3]] <- ts[seq_len(D_cen)]
    dimnames(out_mat_y_pred)[[3]] <- ts[seq_len(D_cen)]

    out_list <- list(mu = out_mat_mu,
                     y_pred = out_mat_y_pred,
                     index = index_df)
    return(out_list)
  }
}
