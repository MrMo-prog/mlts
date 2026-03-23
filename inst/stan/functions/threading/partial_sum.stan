real partial_sum_log_lik(array[] int slice_N, int start, int end,
    array[] int N_obs_id,
    array[] int g_id,
    array[] vector b_free,
    array[] row_vector gammas,
    array[] matrix SIGMA,
    int D_cen,
    int maxLag,
    int D,
    array[] int is_wcen,
    array[] vector y,
    array[] int pos_start,
    array[] int pos_end,
    array[] int pos_start_cov,
    array[] int pos_end_cov,
    matrix b,
    array[] int D_cen_pos,
    array[] int N_pred,
    array[,] int Lag_pred,
    array[,] int D_pred,
    array[,] int D_pred2,
    array[,] int Lag_pred2,
    array[] int Dpos1,
    array[] int Dpos2,
    array[] vector sd_noise,
    // missings and censoring
    int n_miss,
    array[] int n_miss_D,
    array[,] int pos_miss_D,
    vector y_impute,
    array[,] int pos_start_miss,
    array[,] int pos_end_miss,
    array[,] int seq_N_miss,

    int n_censL,
    array[] int n_censL_D,
    array[,] int pos_censL_D,
    vector y_impute_censL,
    array[,] int pos_start_censL,
    array[,] int pos_end_censL,
    array[,] int seq_N_censL,

    int n_censR,
    array[] int n_censR_D,
    array[,] int pos_censR_D,
    vector y_impute_censR,
    array[,] int pos_start_censR,
    array[,] int pos_end_censR,
    array[,] int seq_N_censR,

    int n_inno_covs,
    array[] vector eta_cov,
    array[] int inno_cov_load,
    matrix bmu,
    array[] vector sd_R,
    array[] vector sd_inncov,
    int n_random
    ) {

      real lp = 0;

      for(i in 1:size(slice_N)){
        int pp = slice_N[i];

        int obs_id = N_obs_id[pp]; // observations per person
        int gg_p = g_id[pp];       // group

        // level 2 prediction
        if (n_random == 1){
          lp += normal_lpdf(b_free[pp, 1] | bmu[pp, 1], sd_R[gg_p, 1]);
        } else{
          lp += multi_normal_cholesky_lpdf(b_free[pp, 1: n_random]| bmu[pp, 1:n_random], SIGMA[gg_p]);
        }
        array[n_inno_covs] vector[obs_id - maxLag] eta_cov_id;
        if(n_inno_covs > 0){
          for (n_inno in 1:n_inno_covs){
            eta_cov_id[n_inno, ] = eta_cov[n_inno, pos_start_cov[pp]:pos_end_cov[pp]];
            lp += normal_lpdf(eta_cov_id[n_inno,] | 0, sd_inncov[n_inno, pp]);
          }
        }

        // array of predicted values
        array[D_cen] vector[obs_id-maxLag] mus;
        // create latent mean centered versions of observations
        array[D] vector[N_obs_id[pp]] y_cen;

        for(d in 1:D){
          // calculating missings
          vector[N_obs_id[pp]] y_lokal = y[d, pos_start[pp] : pos_end[pp]]; //slicing person relevant data

          if (seq_N_miss[pp, d] > 0){
            array[seq_N_miss[pp, d]] int y_lokal_pos = pos_miss_D[d, pos_start_miss[pp, d] : pos_end_miss[pp, d]];
            for (m in 1:seq_N_miss[pp, d]){
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute, pos_start_miss[pp, d], seq_N_miss[pp, d]);
          }
          // censoring
          if (seq_N_censL[pp, d] > 0){
            array[seq_N_censL[pp, d]] int y_lokal_pos = pos_censL_D[d, pos_start_censL[pp, d] : pos_end_censL[pp, d]];
            for (m in 1:seq_N_censL[pp, d]){
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute_censL, pos_start_censL[pp, d], seq_N_censL[pp, d]);
          }
          if (seq_N_censR[pp, d] > 0){
            array[seq_N_censR[pp, d]] int y_lokal_pos = pos_censR_D[d, pos_start_censR[pp, d] : pos_end_censR[pp, d]];
            for (m in 1:seq_N_censR[pp, d]){
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute_censR, pos_start_censR[pp, d], seq_N_censR[pp, d]);
          }

          // calculating y_cen --> array vector of within centered observations
          if(is_wcen[d] == 1){
            y_cen[d] = y_lokal - b[pp, D_cen_pos[d]];
          } else {
            y_cen[d] = y_lokal;
          }
        }

        for(d in 1:D){

          if(is_wcen[d] == 1){
            // build prediction matrix for specific dimensions
            int n_cols; // matrix dimensions
            n_cols = (n_inno_covs>0 && d<3) ? N_pred[d] + n_inno_covs : N_pred[d];
            {
              matrix[obs_id - maxLag, n_cols] b_mat; // dimension specific prediction matrix: time points * predictors
              vector[n_cols] b_use; //
              for(nd in 1:N_pred[d]){ // AR effect and CL effects
                int lag_use = Lag_pred[d, nd];
                if(D_pred2[d, nd] == -99){
                  b_mat[,nd] = y_cen[D_pred[d, nd], (1+maxLag-lag_use):(obs_id-lag_use)];
                } else { // interactions between two ds
                  int lag_use2 = Lag_pred2[d, nd];
                  b_mat[,nd] = y_cen[D_pred[d, nd], (1+maxLag-lag_use):(obs_id-lag_use)] .*
                  y_cen[D_pred2[d, nd], (1+maxLag-lag_use2):(obs_id-lag_use2)];
                }
              }
              b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);

              if(n_inno_covs>0 && d < 3){  // add latent factor scores
                b_mat[,(N_pred[d]+1)] = eta_cov_id[1,]; // add innovation covariance factor scores
                b_use[N_pred[d]+1] = inno_cov_load[d];
              }

              mus[D_cen_pos[d]] = b_mat * b_use;
            }
            lp += normal_lpdf(y_cen[d, (1+maxLag):obs_id] | mus[D_cen_pos[d]], sd_noise[D_cen_pos[d],pp]);
          }
        } // end of loop over dimensions

    } // end of loop over subjects
  return(lp);

}
