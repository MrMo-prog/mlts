mlts_residuals <- function(fit_obj, ts, summary = TRUE, se = TRUE, ci = TRUE, ci_level = 0.95) {
  if ("Measurement" %in% fit_obj$model$Model) {
    stop("This function is currently not supported for latent multi-indicator models.", call. = FALSE)
  }
  if (!summary && (se || ci)) {
    warning("Arguments 'ci', 'se' and 'ci_level' will be ignored when summary = FALSE.", call. = FALSE)
  }

  std <- fit_obj$standata
  D_cen <- std$D_cen

  predictions <- mlts_predict(fit_obj = fit_obj, ts = ts, summary = FALSE, ci = FALSE, se = FALSE)
  resid <- predictions$y_cen - predictions$mu
  dimnames(resid)[[3]] <- ts[seq_len(D_cen)]

  if (!summary) {
    return(resid)
  } else {
    
    resid_means <- apply(resid, c(2, 3), mean)
    colnames(resid_means) <- paste0(ts[seq_len(D_cen)], "_resid")
    result_df <- as.data.frame(resid_means)
    if (se) {
      resid_se <- apply(resid, c(2, 3), sd)
      colnames(resid_se) <- paste0(ts[seq_len(D_cen)], "_se")
      result_df <- cbind(result_df, as.data.frame(resid_se))
    }
    if (ci) {
      probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
      q_pp <- apply(resid, c(2, 3), quantile, probs = probs, na.rm = TRUE)
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

# Models fitted with tinterval + missing values, where create_missings() inserts structural NA rows → VARprepare converts these to -Inf in fit_obj$data
# da ist iwas mit den missings falsch idk muss ich ma schauen, anscheinend werden die zu NA und nicht zu -99 oder imputiert oder whatever 
# oder iwie unendlich oder so 
# censored values kriege ich nicht aus standata, da muss ich ma schauen woher ich die bekommen kann oder ob ich die brauche idk

