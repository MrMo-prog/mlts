real partial_sum_log_lik(array[] int slice_N, int start, int end,
    array[] int N_obs_id,
    array[] int g_id,
    array[] vector b_free,
    array[] row_vector gammas,
    array[] matrix SIGMA,
    array[] matrix SIGMA2,
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

    int n_inno_covs,
    array[] vector eta_cov,
    array[] int inno_cov_load,
    matrix bmu,
    array[] vector sd_R,
    array[] vector sd_inncov,
    int n_random,

    //new
    int n_mvn1,
    int n_mvn2,
    int n_iid,
    array[] int pos_iid,
    array[] int pos_mvn1,
    array[] int pos_mvn2,

    // rdsem 
    array[] int is_rdsem,
    array[] int N_pred_rdsem,
    array[,] int D_pred_rdsem,
    array[] int Dpos1_rdsem
    ) {

      real lp = 0;

      for(i in 1:size(slice_N)){
        int pp = slice_N[i];

        int obs_id = N_obs_id[pp]; // observations per person
        int gg_p = g_id[pp];       // group

        // level 2 prediction
        // individual parameters from (multivariate) normal distribution
        if(n_iid > 0){
          for(jj in pos_iid){
            lp += normal_lpdf(b_free[pp,jj] | bmu[pp,jj], sd_R[gg_p, jj]);
          }
        }
        if(n_mvn1 > 0) {
          lp += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn1] | bmu[pp, pos_mvn1], SIGMA[gg_p]);
        }
        if(n_mvn2 > 0) {
          lp += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn2] | bmu[pp, pos_mvn2], SIGMA2[gg_p]);
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
          vector[N_obs_id[pp]] y_lokal = y[d, pos_start[pp] : pos_end[pp]];
          // calculating y_cen --> array vector of within centered observations
          if(is_wcen[d] == 1){
            y_cen[d] = y_lokal - b[pp, D_cen_pos[d]];
          } else {
            y_cen[d] = y_lokal;
          }
        }

        if( sum(is_rdsem) > 0 ){
          for (d in 1:D) {
            if (is_rdsem[d] == 1 ) {
              for(k in 1:N_pred_rdsem[d]){
                y_cen[d,] = y_cen[d,] - b[pp,(Dpos1_rdsem[d]+(k-1))] * y_cen[D_pred_rdsem[d,k],];
              }
            }
          }
        }

        for(d in 1:D){

          if(is_wcen[d] == 1){
            // build prediction matrix for specific dimensions
            int n_cols; // matrix dimensions
            n_cols = (n_inno_covs>0 && d<3) ? N_pred[d] + n_inno_covs : N_pred[d];
            if (N_pred[d] > 0) {
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
            } else {
              mus[D_cen_pos[d],] = rep_vector(0, (obs_id - maxLag));
            }
            lp += normal_lpdf(y_cen[d, (1+maxLag):obs_id] | mus[D_cen_pos[d]], sd_noise[D_cen_pos[d],pp]);
          }
        } // end of loop over dimensions

    } // end of loop over subjects
  return(lp);

}
