predict_iteration <- function(fit_obj, ts, b) {

  D <- fit_obj$standata$D
  N <- fit_obj$standata$N
  N_obs_id <- fit_obj$standata$N_obs_id
  is_wcen <- fit_obj$standata$is_wcen

  N_pred <- fit_obj$standata$N_pred
  maxLag <- fit_obj$standata$maxLag
  Lag_pred <- fit_obj$standata$Lag_pred
  Lag_pred2 <- fit_obj$standata$Lag_pred2
  D_cen <- fit_obj$standata$D_cen
  D_cen_pos <- fit_obj$standata$D_cen_pos
  D_pred2 <- fit_obj$standata$D_pred2
  D_pred <- fit_obj$standata$D_pred
  Dpos1 <- fit_obj$standata$Dpos1
  Dpos2 <- fit_obj$standata$Dpos2

  is_rdsem <- fit_obj$standata$is_rdsem
  N_pred_rdsem <- fit_obj$standata$N_pred_rdsem


  # step 1 --> calcualte y_cen
  y_cen <- vector("list", N)
  pos <- 1

  for (pp in seq_len(N)){
      obs_id = N_obs_id[pp]
      y_cen[[pp]] <- matrix(NA_real_, nrow = obs_id, ncol = D)
    for (d in seq_len(D)){
      if (is_wcen[d] == 1){
        y_cen[[pp]][,d] <- fit_obj$data[pos:(pos + obs_id - 1), ts[d]] - 
          b[pp, D_cen_pos[d]]
      } else {
        y_cen[[pp]][,d] <- fit_obj$data[pos:(pos + obs_id - 1), ts[d]]
      }
    }


    # check for rdsem
    if (sum(is_rdsem) > 0){
      for (d in seq_len(D)){
        if (is_rdsem[d] == 1){
          for (k in seq_len(N_pred_rdsem[d])){
            y_cen[[pp]][, d] <- y_cen[[pp]][, d] - b[pp, (fit_obj$standata$Dpos1_rdsem[d] + (k-1))] * y_cen[[pp]][, fit_obj$standata$D_pred_rdsem[d, k]]
          }
        }
      }
    }
    pos <- pos + obs_id
  }

  # step 2 --> calculate mus
  mus <- vector("list", N)
  for (pp in seq_len(N)){
    obs_id = N_obs_id[pp]
    mus[[pp]] <- matrix(NA_real_, nrow = obs_id - maxLag, ncol = D_cen)
    for (d in seq_len(D)){
      if (is_wcen[d] == 1){
#        n_cols <- ifelse(fit_obj$standata$n_inno_covs > 0 && d < 3,
#        N_pred[d] + fit_obj$standata$n_inno_covs, N_pred[d])
        n_cols <- N_pred[d]

        if (N_pred[d] > 0){
          b_mat <- matrix(NA_real_, nrow = (obs_id - maxLag), ncol = n_cols)
          b_use <- numeric(n_cols)
          for (nd in seq_len(N_pred[d])){
            lag_use <- Lag_pred[d, nd]
            if (D_pred2[d, nd]  == -99){
              b_mat[, nd] <- y_cen[[pp]][(1 + maxLag - lag_use) : (obs_id - lag_use), D_pred[d, nd]]
            } else {
              lag_use2 <- Lag_pred2[d, nd]
              b_mat[, nd] <- y_cen[[pp]][(1 + maxLag - lag_use) : (obs_id - lag_use), D_pred[d, nd]] *
              y_cen[[pp]][(1 + maxLag - lag_use2) : (obs_id - lag_use2), D_pred2[d, nd]]
            }
          }
          b_use[1:N_pred[d]] <- b[pp, Dpos1[d]:Dpos2[d]]

#          if (fit_obj$standata$n_inno_covs > 0 && d < 3){
#            b_mat[,(N_pred[d] + 1 )] <- fit_obj$standata$eta_cov_id[1, ]
#            b_use[N_pred[d] + 1] <- fit_obj$standata$inno_cov_load[d]
#          }

          mus[[pp]][, D_cen_pos[d]] <- b_mat %*% b_use
        } else {
          mus[[pp]][, D_cen_pos[d]] <- rep(0, obs_id - maxLag)
        }
      }
    }
  }

  return(list(mus = mus, y_cen = y_cen))
}