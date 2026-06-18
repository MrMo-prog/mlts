real partial_sum_log_lik(array[] vector slice_b_free, int start, int end,
    array[] int N_obs_id,
    array[] int g_id,
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
    array[] int D_cen_pos,
    array[] int N_pred,
    array[,] int Lag_pred,
    array[,] int D_pred,
    array[,] int D_pred2,
    array[,] int Lag_pred2,
    array[] int Dpos1,
    array[] int Dpos2,
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
    array[] vector sd_R,
    int n_random,

    // new
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
    array[] int Dpos1_rdsem,
    
    // parameters to calculate b, bmu, sd_noise, sd_inncov locally
    int n_pars,
    array[] int is_random,
    int n_fixed,
    array[] vector b_fix,
    array[,] int is_fixed,
    int n_cov,
    int n_cov_bs,
    array[,] int n_cov_mat,
    array[] vector b_re_pred,
    matrix W,
    array[] int innos_rand,
    array[] int innos_fix_pos,
    array[] int innos_pos,
    array[,] int inno_cov_pos,
    array[] vector sigma
    ) {

      real lp = 0;

      for (i in 1:size(slice_b_free)) {
        int pp = start + i - 1; // absolute subject index

        int obs_id = N_obs_id[pp]; // observations per person
        int gg_p = g_id[pp];       // group
        vector[n_random] b_free_pp = slice_b_free[i];

        // 1. Calculate bmu_pp (size n_random)
        row_vector[n_random] bmu_pp;
        {
          matrix[n_cov, n_random] b_re_pred_mat_g = rep_matrix(0, n_cov, n_random);
          b_re_pred_mat_g[1, ] = gammas[gg_p, ];
          if (n_cov > 1) {
            for (idx in 1:n_cov_bs) {
              b_re_pred_mat_g[n_cov_mat[idx, 1], n_cov_mat[idx, 2]] = b_re_pred[gg_p, idx];
            }
          }
          bmu_pp = W[pp, ] * b_re_pred_mat_g;
        }

        // level 2 prediction
        // individual parameters from (multivariate) normal distribution
        if (n_iid > 0) {
          for (jj in pos_iid) {
            lp += normal_lpdf(b_free_pp[jj] | bmu_pp[jj], sd_R[gg_p, jj]);
          }
        }
        if (n_mvn1 > 0) {
          lp += multi_normal_cholesky_lpdf(b_free_pp[pos_mvn1] | to_vector(bmu_pp[pos_mvn1]), SIGMA[gg_p]);
        }
        if (n_mvn2 > 0) {
          lp += multi_normal_cholesky_lpdf(b_free_pp[pos_mvn2] | to_vector(bmu_pp[pos_mvn2]), SIGMA2[gg_p]);
        }

        // 2. Calculate b_pp (size n_pars)
        vector[n_pars] b_pp;
        for (j in 1:n_random) {
          b_pp[is_random[j]] = b_free_pp[j];
        }
        if (n_fixed > 0) {
          for (j in 1:n_fixed) {
            b_pp[is_fixed[1, j]] = b_fix[gg_p, j];
          }
        }

        // 3. Calculate sd_inncov_pp (size n_inno_covs)
        vector[n_inno_covs] sd_inncov_pp;
        if (n_inno_covs > 0) {
          for (idx in 1:n_inno_covs) {
            sd_inncov_pp[idx] = sqrt(exp(b_pp[inno_cov_pos[1, idx]]));
          }
        }

        array[n_inno_covs] vector[obs_id - maxLag] eta_cov_id;
        if (n_inno_covs > 0) {
          for (n_inno in 1:n_inno_covs) {
            eta_cov_id[n_inno, ] = eta_cov[n_inno, pos_start_cov[pp]:pos_end_cov[pp]];
            lp += normal_lpdf(eta_cov_id[n_inno, ] | 0, sd_inncov_pp[n_inno]);
          }
        }

        // array of predicted values
        array[D_cen] vector[obs_id - maxLag] mus;
        // create latent mean centered versions of observations
        array[D] vector[N_obs_id[pp]] y_cen;

        for (d in 1:D) {
          int offset_miss = 0;
          int offset_censL = 0;
          int offset_censR = 0;
          if (d > 1) {
            offset_miss = sum(n_miss_D[1 : (d - 1)]);
            offset_censL = sum(n_censL_D[1 : (d - 1)]);
            offset_censR = sum(n_censR_D[1 : (d - 1)]);
          }
          // calculating missings
          vector[N_obs_id[pp]] y_lokal = y[d, pos_start[pp] : pos_end[pp]]; // slicing person relevant data

          if (seq_N_miss[pp, d] > 0) {
            array[seq_N_miss[pp, d]] int y_lokal_pos = pos_miss_D[d, (pos_start_miss[pp, d] - offset_miss) : (pos_end_miss[pp, d] - offset_miss)];
            for (m in 1:seq_N_miss[pp, d]) {
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute, pos_start_miss[pp, d], seq_N_miss[pp, d]);
          }
          // censoring
          if (seq_N_censL[pp, d] > 0) {
            array[seq_N_censL[pp, d]] int y_lokal_pos = pos_censL_D[d, (pos_start_censL[pp, d] - offset_censL) : (pos_end_censL[pp, d] - offset_censL)];
            for (m in 1:seq_N_censL[pp, d]) {
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute_censL, pos_start_censL[pp, d], seq_N_censL[pp, d]);
          }
          if (seq_N_censR[pp, d] > 0) {
            array[seq_N_censR[pp, d]] int y_lokal_pos = pos_censR_D[d, (pos_start_censR[pp, d] - offset_censR) : (pos_end_censR[pp, d] - offset_censR)];
            for (m in 1:seq_N_censR[pp, d]) {
              y_lokal_pos[m] = y_lokal_pos[m] - pos_start[pp] + 1;
            }
            y_lokal[y_lokal_pos] = segment(y_impute_censR, pos_start_censR[pp, d], seq_N_censR[pp, d]);
          }

          // calculating y_cen --> array vector of within centered observations
          if (is_wcen[d] == 1) {
            y_cen[d] = y_lokal - b_pp[D_cen_pos[d]];
          } else {
            y_cen[d] = y_lokal;
          }
        }

        if (sum(is_rdsem) > 0) {
          for (d in 1:D) {
            if (is_rdsem[d] == 1) {
              for (k in 1:N_pred_rdsem[d]) {
                y_cen[d] = y_cen[d] - b_pp[(Dpos1_rdsem[d] + (k - 1))] * y_cen[D_pred_rdsem[d, k]];
              }
            }
          }
        }

        // 4. Calculate sd_noise_pp (size D_cen)
        vector[D_cen] sd_noise_pp;
        for (d in 1:D_cen) {
          if (innos_rand[d] == 0) {
            sd_noise_pp[d] = sigma[gg_p, innos_fix_pos[d]];
          } else {
            sd_noise_pp[d] = sqrt(exp(b_pp[innos_pos[d]]));
          }
        }

        for (d in 1:D) {
          if (is_wcen[d] == 1) {
            // build prediction matrix for specific dimensions
            int n_cols; // matrix dimensions
            n_cols = (n_inno_covs > 0 && d < 3) ? N_pred[d] + n_inno_covs : N_pred[d];
            if (N_pred[d] > 0) {
              matrix[obs_id - maxLag, n_cols] b_mat; // dimension specific prediction matrix: time points * predictors
              vector[n_cols] b_use;
              for (nd in 1:N_pred[d]) { // AR effect and CL effects
                int lag_use = Lag_pred[d, nd];
                if (D_pred2[d, nd] == -99) {
                  b_mat[ , nd] = y_cen[D_pred[d, nd], (1 + maxLag - lag_use):(obs_id - lag_use)];
                } else { // interactions between two ds
                  int lag_use2 = Lag_pred2[d, nd];
                  b_mat[ , nd] = y_cen[D_pred[d, nd], (1 + maxLag - lag_use):(obs_id - lag_use)] .*
                                 y_cen[D_pred2[d, nd], (1 + maxLag - lag_use2):(obs_id - lag_use2)];
                }
              }
              b_use[1:N_pred[d]] = b_pp[Dpos1[d]:Dpos2[d]];

              if (n_inno_covs > 0 && d < 3) { // add latent factor scores
                b_mat[ , (N_pred[d] + 1)] = eta_cov_id[1]; // add innovation covariance factor scores
                b_use[N_pred[d] + 1] = inno_cov_load[d];
              }

              mus[D_cen_pos[d]] = b_mat * b_use;
            } else {
              mus[D_cen_pos[d]] = rep_vector(0, (obs_id - maxLag));
            }
            lp += normal_lpdf(y_cen[d, (1 + maxLag):obs_id] | mus[D_cen_pos[d]], sd_noise_pp[D_cen_pos[d]]);
          }
        } // end of loop over dimensions

      } // end of loop over subjects
      return lp;
}
