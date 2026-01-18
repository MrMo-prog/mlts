array[] vector calculate_mus(
  array[] vector x_dyn, // y_cen oder etaW_id
  int is_latent, // 1=latent, 0=manifest,
  int is_covs_fix, // 0=nein, 1=covs_fix

  // indices / sizes
  int pp,
  int obs_id,
  int maxLag,
  int D,
  int D_cen,

  // mapping & model structure
  array[] int is_wcen,
  array[] int D_cen_pos,
  array[] int N_pred,
  array[,] int D_pred,
  array[,] int Lag_pred,
  array[,] int D_pred2,
  array[,] int Lag_pred2,
  array[] int Dpos1,
  array[] int Dpos2,

  // parameters
  matrix b,

  // innovation cov --> dummy code for covsfix models
  int n_inno_covs,
  array[] int inno_cov_load,
  array[] vector eta_cov_id,

  // latent SI intercept helpers --> dummy code for manifest models and latent covsfix
  array[] int D_np,
  array[] int D_pos_is_SI,
  array[] vector YB){

  int T = obs_id - maxLag;
  matrix[T, D_cen] mus;

  for(d in 1:D){ // start loop over dimensions

      if(is_wcen[d] == 1){

        // build prediction matrix for specific dimensions
        int n_cols = N_pred[d];  // matrix dimensions
        matrix[(T),n_cols] b_mat;
        vector[n_cols] b_use;

      for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
         int lag_use = Lag_pred[d,nd];
         if(D_pred2[d,nd] == -99){
          b_mat[,nd] = x_dyn[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)];
         } else {
          int lag_use2 = Lag_pred2[d,nd];
          b_mat[,nd] = x_dyn[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)] .*
                       x_dyn[D_pred2[d, nd],(1+maxLag-lag_use2):(obs_id-lag_use2)];
         }
      }
      b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);
      vector[T] b_calculated;
      // manifest
      if (is_latent == 0){
        b_calculated = b[pp,D_cen_pos[d]] + b_mat * b_use;
      }
      // latent
      else if (is_latent == 1){
        b_calculated = b_mat * b_use;
      }
      // innocov dazu addieren, wenn =! covsfix
      if (is_covs_fix == 0 && n_inno_covs > 0 && d < 3){ // !!!hardcode --> nur 1 cov erlaubt
        b_calculated += inno_cov_load[d] * eta_cov_id[1]; // !!!hardcode --> nur 1 cov erlaubt
      }
      //
      if ((D_np[d] == 1) && is_latent == 1 && is_covs_fix == 0){
          b_calculated += YB[D_pos_is_SI[d],pp];
        }
      mus[, D_cen_pos[d]] = b_calculated;
    }
  }
  if (is_covs_fix == 1){
    array[T] vector[D_cen] out;
    for (t in 1:T) out[t] = mus[t]';
    return out;

  } else {
      array[D_cen] vector[T] out;
      for (d in 1:D_cen) out[d] = col(mus, d);
      return out;
  }
}

