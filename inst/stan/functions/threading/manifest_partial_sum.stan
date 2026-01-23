real partial_sum_manifest (array[] int seq_slice, int start, int end,
// Größen & Indizes
    int N,
    array[] int N_obs_id,
    array[] int g_id,
    int maxLag,
    int D,
    int D_cen,
    int n_random,
    int n_inno_covs,

    // 2. Daten & Parameter
    array[] vector y_merge,      // Die Beobachtungen
    matrix bmu,                  // Transformed Parameter
    array[] vector sd_R,         // Parameter
    array[] matrix SIGMA,        // ACHTUNG: Das wird im Model-Block berechnet!
    array[] vector b_free,       // Parameter (für das Sampling)
    matrix b,                    // Transformed Parameter
    array[] vector eta_cov,      // Parameter
    array[] vector sd_inncov,    // Transformed Parameter
    array[] vector sd_noise,     // Transformed Parameter

    // --- 3. Argumente für calculate_mus ---
    array[] int is_wcen,
    array[] int D_cen_pos,
    array[] int N_pred,
    array[,] int D_pred,         // 2D Array
    array[,] int Lag_pred,       // 2D Array
    array[,] int D_pred2,        // 2D Array
    array[,] int Lag_pred2,      // 2D Array
    array[] int Dpos1,
    array[] int Dpos2,
    array[] int inno_cov_load,

    // Dummys
    array[] int dummy_D_np,
    array[] int dummy_D_pos_is_SI,
    array[] vector dummy_YB,

    // 4. Index-Arrays (fürs Threading)
    array[] int pos_start,
    array[] int pos_cov_start){

  real lp = 0;

  for (j in 1:size(seq_slice)) {

    int pp = seq_slice[j];
    // store number of observations per person
    int obs_id = (N_obs_id[pp]);
    // obtain group
    int pp_g = g_id[pp];

    int pos = pos_start[pp];
    int pos_cov = pos_cov_start[pp];

    // individual parameters from (multivariate) normal distribution
    if(n_random == 1){
      lp += normal_lpdf(b_free[pp,1] | bmu[pp,1], sd_R[pp_g, 1]);
    } else {
      lp += multi_normal_cholesky_lpdf(b_free[pp, 1:n_random] | bmu[pp, 1 : n_random], SIGMA[pp_g]);
    }

    array[n_inno_covs] vector[obs_id-maxLag] eta_cov_id;
    if(n_inno_covs>0){
      for(i in 1:n_inno_covs){
          eta_cov_id[i,] = segment(eta_cov[i,], pos_cov, (obs_id-maxLag));
          lp += normal_lpdf(eta_cov_id[i,] | 0, sd_inncov[i,pp]);
      }
    }

    // local variable declaration: array of predicted values
    {
    array[D_cen] vector[obs_id-maxLag] mus;

    // create latent mean centered versions of observations
    array[D] vector[obs_id] y_cen;

    for(d in 1:D){ // start loop over dimensions
      if(is_wcen[d] == 1){
        y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)] - b[pp,D_cen_pos[d]];
      } else {
        y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)];
      }
    }

    mus = calculate_mus(
            y_cen,              // x_dyn (here: y_cen)
            0,                  // is_latent = 0 (Manifest)
            0,                  // is_covs_fix = 0 (Manifest Covs)
            pp, obs_id, maxLag, D, D_cen, is_wcen, D_cen_pos, N_pred, D_pred,
            Lag_pred, D_pred2, Lag_pred2, Dpos1, Dpos2, b, n_inno_covs,
            inno_cov_load, eta_cov_id, dummy_D_np, dummy_D_pos_is_SI, dummy_YB);

    for(d in 1:D){ // start loop over dimensions
      if(is_wcen[d] == 1){
        // sampling statement
        lp += normal_lpdf(y_merge[d,(pos+maxLag):(pos+(obs_id-1))] | mus[D_cen_pos[d],], sd_noise[D_cen_pos[d],pp]);
       }
      }
    }
  } // end loop over subjects
  return lp;
}


