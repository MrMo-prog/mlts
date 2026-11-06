mlts_residuals <- function(fit_obj, ts, summary = TRUE, se = TRUE, ci = TRUE, ci_level = 0.95) {
  if ("Measurement" %in% fit_obj$model$Model) {
    stop("This function is currently not supported for latent multi-indicator models.", call. = FALSE)
  }
  if (!summary && (se || ci)) {
    warning("Arguments 'ci', 'se' and 'ci_level' will be ignored when summary = FALSE.", call. = FALSE)
  }

  std <- fit_obj$standata
  D_cen <- std$D_cen
  D <- std$D

  y_obs_mat <- matrix(
  NA_real_,
  nrow = nrow(fit_obj$data),
  ncol = D_cen
  ) 
  is_wcen <- std$is_wcen
  D_cen_pos <- std$D_cen_pos

  for (d in seq_len(D)) {
    if (is_wcen[d] == 1) {
      y_obs_mat[, D_cen_pos[d]] <- fit_obj$data[, ts[d]]
    }
  }

  predictions <- mlts_predict(fit_obj = fit_obj, ts = ts, summary = FALSE, ci = FALSE, se = FALSE)
  resid <- array(NA_real_, dim = c(dim(predictions$y_pred)[1], nrow(fit_obj$data), D_cen ))
  for (s in seq_len(dim(predictions$y_pred)[1])) {
    pred_s <- matrix(
      predictions$y_pred[s, , ],
      nrow = nrow(fit_obj$data),
      ncol = D_cen
    )
    resid[s, , ] <- y_obs_mat - pred_s
  }
  dimnames(resid)[[3]] <- ts[seq_len(D_cen)]
  if (!summary) {
    out_list <- list(residuals = resid, index = predictions$index)
    return(out_list)
  } else {

    result_df <- predictions$index
    resid_means <- apply(resid, c(2, 3), mean)
    colnames(resid_means) <- paste0(ts[seq_len(D_cen)], "_resid")
    result_df <- cbind(result_df, as.data.frame(resid_means))
    if (se) {
      resid_se <- apply(resid, c(2, 3), sd)
      colnames(resid_se) <- paste0(ts[seq_len(D_cen)], "_se")
      result_df <- cbind(result_df, as.data.frame(resid_se))
    }
    if (ci) {
      probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
      q_pp <- apply(resid, c(2, 3), function(x) {
        if (all(is.na(x))) {
          return(c(NA_real_, NA_real_))
        }
        quantile(x, probs = probs, na.rm = TRUE)
      })
      dim(q_pp) <- c(2, dim(resid_means))

      resid_lower <- array(q_pp[1, , , drop = FALSE], dim = dim(resid_means))
      resid_upper <- array(q_pp[2, , , drop = FALSE], dim = dim(resid_means))
      colnames(resid_lower) <- paste0(ts[seq_len(D_cen)], "_lower")
      colnames(resid_upper) <- paste0(ts[seq_len(D_cen)], "_upper")
      result_df <- cbind(result_df, as.data.frame(resid_lower), as.data.frame(resid_upper))
    }
    return(result_df)

  }
}
