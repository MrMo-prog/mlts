mlts_predict <- function(fit_obj, ts, summary = TRUE, se = TRUE, ci = TRUE, ci_level = 0.95){
  std <- fit_obj$standata

  if ("Measurement" %in% fit_obj$model$Model) {
    stop("This function is currently not supported for latent multi-indicator models.", call. = FALSE)
  }
  if (!summary && (se || ci)) {
    warning("Arguments 'ci', 'se' and 'ci_level' will be ignored when summary = FALSE.", call. = FALSE)
  }
  has_b_free <- any(grepl("^b_free", names(fit_obj$stanfit)))
  if (std$n_random > 0 && !has_b_free) {
    stop("This method currently requires 'monitor_person_pars = TRUE' during fitting.")
  }
  # TODO: Implement fallback for monitor_person_pars = FALSE using person.pars.summary


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
  y_cen_all <- vector("list", N)

  for (pp in seq_len(N)){
    obs_p <- N_obs_id[pp]
    mus_all[[pp]] <- array(NA_real_, dim = c(S, obs_p - maxLag, D_cen))
    y_cen_all[[pp]] <- array(NA_real_, dim = c(S, obs_p, D))
  }
  for (s in seq_len(S)){
    result <- predict_iteration(fit_obj, ts, b[s, , ])

    for (pp in seq_len(N)){
      mus_all[[pp]][s, , ] <- result$mus[[pp]]
      y_cen_all[[pp]][s, , ] <- result$y_cen[[pp]]
    }
  }
  if (summary == TRUE){
    mu_means <- vector("list", N)
    if (se) {
      mu_ses <- vector("list", N)
      pred_mat_se <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    }
    if (ci) {
      mu_lower <- vector("list", N)
      mu_upper <- vector("list", N)
      pred_mat_lower <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
      pred_mat_upper <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    }
    if (S > 1){
      probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
      for (pp in seq_len(N)){
        mu_means[[pp]] <- apply(mus_all[[pp]], c(2, 3), mean)
        if (se == TRUE){
          mu_ses[[pp]] <- apply(mus_all[[pp]], c(2, 3), sd)
        }
        if (ci == TRUE) {
          q_pp <- apply(mus_all[[pp]], c(2, 3), quantile, probs = probs, na.rm = TRUE)
          mu_lower[[pp]] <- array(q_pp[1, , ], dim = dim(mus_all[[pp]])[2:3])
          mu_upper[[pp]] <- array(q_pp[2, , ], dim = dim(mus_all[[pp]])[2:3])
        }
      }
    } else {
      for (pp in seq_len(N)) {
        mu_means[[pp]] <- array(mus_all[[pp]][1, , , drop = FALSE], dim = dim(mus_all[[pp]])[2:3])
      }
    }
    pred_mat_mu <- matrix(NA_real_, nrow = nrow(fit_obj$data), ncol = D_cen)
    pos = 1
    for (pp in seq_len(N)){
      n_obs <- N_obs_id[pp]
      pred_mat_mu[(pos + maxLag):(pos + n_obs - 1), ] <- mu_means[[pp]]
      if (se) pred_mat_se[(pos + maxLag):(pos + n_obs - 1), ] <- mu_ses[[pp]]
      if (ci) {
        pred_mat_lower[(pos + maxLag):(pos + n_obs - 1), ] <- mu_lower[[pp]]
        pred_mat_upper[(pos + maxLag):(pos + n_obs - 1), ] <- mu_upper[[pp]]
      }
      pos <- pos + n_obs
    }
    result_df <- as.data.frame(pred_mat_mu)
    colnames(result_df) <- paste0(ts[seq_len(D_cen)], "_pred")

    if (se) {
      se_df <- as.data.frame(pred_mat_se)
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

  } else { 
    out_mat_mu <- array(NA_real_, dim = c(S, nrow(fit_obj$data), D_cen))
    out_mat_y_cen <- array(NA_real_, dim = c(S, nrow(fit_obj$data), D_cen))
    pos = 1
    for (pp in seq_len(N)){
      n_obs <- N_obs_id[pp]
      out_mat_mu[, (pos + maxLag):(pos + n_obs - 1), ] <- mus_all[[pp]]
      out_mat_y_cen[, pos:(pos + n_obs - 1), ] <- y_cen_all[[pp]][, , seq_len(D_cen)]
      pos = pos + n_obs
    }
    dimnames(out_mat_mu)[[3]] <- ts[seq_len(D_cen)]
    dimnames(out_mat_y_cen)[[3]] <- ts[seq_len(D_cen)]
    return(list(mu = out_mat_mu, y_cen = out_mat_y_cen))
  }
}
