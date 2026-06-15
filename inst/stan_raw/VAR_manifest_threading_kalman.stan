functions{
  #include "../stan/functions/function_calculate_b.stan"
  #include "../stan/functions/function_hyper_priors.stan"
  #include "../stan/functions/function_calculate_bmu.stan"
  #include "../stan/functions/function_outcome_prediction.stan"

  real partial_sum_log_lik_kalman(array[] int seq_N_slice,
                           int start, int end,
                           array[] int N_obs_id, array[] int g_id, array[] vector b_free,
                           array[] row_vector gammas, array[] matrix SIGMA, array[] matrix SIGMA2,
                           int D_cen, int maxLag, int D, array[] int is_wcen,
                           array[] vector y, array[] int pos_start, array[] int pos_end,
                           array[] int pos_start_cov, array[] int pos_end_cov, matrix b,
                           array[] int D_cen_pos, array[] int N_pred, array[,] int Lag_pred,
                           array[,] int D_pred, array[,] int D_pred2, array[,] int Lag_pred2,
                           array[] int Dpos1, array[] int Dpos2, array[] vector sd_noise,
                           int n_miss, array[] int n_miss_D, array[,] int pos_miss_D,
                           array[,] int pos_start_miss, array[,] int pos_end_miss, array[,] int seq_N_miss,
                           int n_censL, array[] int n_censL_D, array[,] int pos_censL_D,
                           array[,] int pos_start_censL, array[,] int pos_end_censL, array[,] int seq_N_censL,
                           int n_censR, array[] int n_censR_D, array[,] int pos_censR_D,
                           array[,] int pos_start_censR, array[,] int pos_end_censR, array[,] int seq_N_censR,
                           int n_inno_covs, array[] vector eta_cov, array[] int inno_cov_load,
                           matrix bmu, array[] vector sd_R, array[] vector sd_inncov,
                           int n_random, int n_mvn1, int n_mvn2, int n_iid, array[] int pos_iid,
                           array[] int pos_mvn1, array[] int pos_mvn2,
                           array[] int is_rdsem, array[] int N_pred_rdsem, array[,] int D_pred_rdsem,
                           array[] int Dpos1_rdsem, real censL_val, real censR_val) {
    real pt = 0;
    
    for(n in 1:size(seq_N_slice)) {
      int pp = seq_N_slice[n];
      int obs_id = N_obs_id[pp];
      int pp_g = g_id[pp];

      if(n_iid > 0){
        for(jj in pos_iid){
          pt += normal_lpdf(b_free[pp,jj] | bmu[pp,jj], sd_R[pp_g, jj]);
        }
      }
      if(n_mvn1 > 0) {
        pt += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn1] | bmu[pp, pos_mvn1], SIGMA[pp_g]);
      }
      if(n_mvn2 > 0) {
        pt += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn2] | bmu[pp, pos_mvn2], SIGMA2[pp_g]);
      }

      array[n_inno_covs] vector[obs_id-maxLag] eta_cov_id;
      if (n_inno_covs > 0){
        for (i in 1:n_inno_covs){
          eta_cov_id[i, ] = eta_cov[i, pos_start_cov[pp] : pos_end_cov[pp]];
          pt += normal_lpdf(eta_cov_id[i, ] | 0, sd_inncov[i, pp]);
        }
      }

      {

      // Kalman Filter setup
      int N_X = D_cen * maxLag;
      vector[N_X] x_upd = rep_vector(0, N_X);
      matrix[N_X, N_X] P_upd = rep_matrix(0, N_X, N_X);

      // Boolean arrays for missing and censoring
      array[D, obs_id] int is_miss = rep_array(0, D, obs_id);
      array[D, obs_id] int is_cL = rep_array(0, D, obs_id);
      array[D, obs_id] int is_cR = rep_array(0, D, obs_id);
      
      for(d in 1:D) {
         for(i in 1:n_miss_D[d]) {
            int g_idx = pos_miss_D[d, i];
            if (g_idx >= pos_start[pp] && g_idx <= pos_end[pp]) {
               is_miss[d, g_idx - pos_start[pp] + 1] = 1;
            }
         }
         for(i in 1:n_censL_D[d]) {
            int g_idx = pos_censL_D[d, i];
            if (g_idx >= pos_start[pp] && g_idx <= pos_end[pp]) {
               is_cL[d, g_idx - pos_start[pp] + 1] = 1;
            }
         }
         for(i in 1:n_censR_D[d]) {
            int g_idx = pos_censR_D[d, i];
            if (g_idx >= pos_start[pp] && g_idx <= pos_end[pp]) {
               is_cR[d, g_idx - pos_start[pp] + 1] = 1;
            }
         }
      }

      // Compute offset (mean centering + RDSEM)
      matrix[D, obs_id] d_offset = rep_matrix(0, D, obs_id);
      for (d in 1:D) {
         if (is_wcen[d] == 1) {
            for (t in 1:obs_id) {
               d_offset[d, t] += b[pp, D_cen_pos[d]];
            }
         }
      }
      if (sum(is_rdsem) > 0) {
         for (d in 1:D) {
            if (is_rdsem[d] == 1) {
               for (k in 1:N_pred_rdsem[d]) {
                  int d_pred = D_pred_rdsem[d, k];
                  real b_val = b[pp, Dpos1_rdsem[d] + k - 1];
                  for (t in 1:obs_id) {
                     // Evaluate at raw exogenous data
                     d_offset[d, t] += b_val * y[d_pred, pos_start[pp] + t - 1];
                  }
               }
            }
         }
      }

      // Observed (adjusted) data matrix
      matrix[D, obs_id] y_cen_obs;
      for (d in 1:D) {
         for (t in 1:obs_id) {
            y_cen_obs[d, t] = y[d, pos_start[pp] + t - 1] - d_offset[d, t];
         }
      }

      // Kalman initialization (first maxLag points perfectly observed)
      for (d in 1:D) {
         if (is_wcen[d] == 1) {
            for (l in 1:maxLag) {
               x_upd[D_cen_pos[d] + (l-1)*D_cen] = y_cen_obs[d, maxLag - l + 1];
            }
         }
      }

      // Kalman sequence
      for (t in (maxLag + 1):obs_id) {
         vector[N_X] x_pred = rep_vector(0, N_X);
         matrix[N_X, N_X] J = rep_matrix(0, N_X, N_X);

         // 1. Shift states
         for (d in 1:D) {
            if (is_wcen[d] == 1) {
               int d_idx = D_cen_pos[d];
               for (l in 2:maxLag) {
                  x_pred[d_idx + (l-1)*D_cen] = x_upd[d_idx + (l-2)*D_cen];
                  J[d_idx + (l-1)*D_cen, d_idx + (l-2)*D_cen] = 1.0;
               }
            }
         }

         // 2. Predict new states
         for (d in 1:D) {
            if (is_wcen[d] == 1) {
               int d_idx = D_cen_pos[d];
               
               if (N_pred[d] > 0) {
                  for (nd in 1:N_pred[d]) {
                     int l1 = Lag_pred[d, nd];
                     int d1 = D_pred[d, nd];
                     real val1;
                     int idx1 = -1;
                     
                     if (is_wcen[d1] == 1) {
                        idx1 = D_cen_pos[d1] + (l1-1)*D_cen;
                        val1 = x_upd[idx1];
                     } else {
                        val1 = y_cen_obs[d1, t - l1];
                     }
                     
                     real b_val = b[pp, Dpos1[d] + nd - 1];
                     
                     if (D_pred2[d, nd] == -99) {
                        x_pred[d_idx] += b_val * val1;
                        if (idx1 != -1) {
                           J[d_idx, idx1] += b_val;
                        }
                     } else {
                        // Interaction
                        int l2 = Lag_pred2[d, nd];
                        int d2 = D_pred2[d, nd];
                        real val2;
                        int idx2 = -1;
                        
                        if (is_wcen[d2] == 1) {
                           idx2 = D_cen_pos[d2] + (l2-1)*D_cen;
                           val2 = x_upd[idx2];
                        } else {
                           val2 = y_cen_obs[d2, t - l2];
                        }
                        
                        x_pred[d_idx] += b_val * val1 * val2;
                        
                        if (idx1 != -1 && idx2 != -1) {
                           x_pred[d_idx] += b_val * P_upd[idx1, idx2];
                           J[d_idx, idx1] += b_val * val2;
                           J[d_idx, idx2] += b_val * val1;
                        } else if (idx1 != -1) {
                           J[d_idx, idx1] += b_val * val2;
                        } else if (idx2 != -1) {
                           J[d_idx, idx2] += b_val * val1;
                        }
                     }
                  }
               }
               
               if (n_inno_covs > 0 && d < 3) {
                  x_pred[d_idx] += inno_cov_load[d] * eta_cov_id[1, t - maxLag];
               }
            }
         }

         // 3. Covariance prediction
         matrix[N_X, N_X] P_pred = J * P_upd * J';
         for (d in 1:D) {
            if (is_wcen[d] == 1) {
               int d_idx = D_cen_pos[d];
               P_pred[d_idx, d_idx] += square(sd_noise[d_idx, pp]);
            }
         }

         // Determine which variables are exactly observed at time t
         int n_obs_exact = 0;
         array[D_cen] int obs_idx;
         array[D_cen] real obs_val;

         for (d in 1:D) {
            if (is_wcen[d] == 1) {
               if (is_miss[d, t] == 0 && is_cL[d, t] == 0 && is_cR[d, t] == 0) {
                  n_obs_exact += 1;
                  obs_idx[n_obs_exact] = D_cen_pos[d];
                  obs_val[n_obs_exact] = y_cen_obs[d, t];
               }
            }
         }

         // Update step
         if (n_obs_exact > 0) {
            vector[n_obs_exact] v_err;
            matrix[n_obs_exact, n_obs_exact] S;
            
            for (i in 1:n_obs_exact) {
               v_err[i] = obs_val[i] - x_pred[obs_idx[i]];
               for (j in 1:n_obs_exact) {
                  S[i, j] = P_pred[obs_idx[i], obs_idx[j]];
               }
            }
            
            for (i in 1:n_obs_exact) {
               S[i, i] += 1e-9;
            }
            
            matrix[n_obs_exact, n_obs_exact] L_S = cholesky_decompose(S);
            pt += multi_normal_cholesky_lpdf(to_vector(obs_val[1:n_obs_exact]) | x_pred[obs_idx[1:n_obs_exact]], L_S);
            
            matrix[N_X, n_obs_exact] PHt;
            for (i in 1:N_X) {
               for (j in 1:n_obs_exact) {
                  PHt[i, j] = P_pred[i, obs_idx[j]];
               }
            }
            
            matrix[N_X, n_obs_exact] K = mdivide_right_spd(PHt, S);
            x_upd = x_pred + K * v_err;
            P_upd = P_pred - K * PHt';
            P_upd = 0.5 * (P_upd + P_upd');
            for (i in 1:N_X) P_upd[i, i] += 1e-9;
         } else {
            x_upd = x_pred;
            P_upd = P_pred;
         }

         // Censored observations contribution
         for (d in 1:D) {
            if (is_wcen[d] == 1) {
               int d_idx = D_cen_pos[d];
               if (is_cL[d, t] == 1) {
                  pt += normal_lcdf(censL_val | x_upd[d_idx] + d_offset[d, t], sqrt(P_upd[d_idx, d_idx]));
               }
               if (is_cR[d, t] == 1) {
                  pt += normal_lccdf(censR_val | x_upd[d_idx] + d_offset[d, t], sqrt(P_upd[d_idx, d_idx]));
               }
            }
         }
      }

      }
    }
    return pt;
  }

}

data {
  int<lower=0> N;                 // number of observational units
  int<lower=1> G;                 // groups in total
  int<lower=1> D; 	              // number of time-varying constructs
  int<lower=1> D_cen;             // number of constructs to be mean-centered
  int<lower=1, upper=3> maxLag;   // maximum lag
  int<lower=1> N_obs; 	          // observations in total: N * TP
  int<lower=1> n_pars;            // number of parameters
  int<lower=1> n_random;          // number of random effects

  //new
  int<lower=0> n_mvn1;
  int<lower=0> n_mvn2;
  int<lower=0> n_iid;
  array[n_iid] int pos_iid;
  array[n_mvn1] int pos_mvn1;
  array[n_mvn2] int pos_mvn2;

  int n_fixed;
  array[1, n_fixed] int is_fixed;
  array[n_random] int is_random;  // which parameters to model person-specific
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[D] vector[N_obs] y; 	    // array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  array[D] int n_miss_D;           // missings per D
  array[D,max(n_miss_D)] int pos_miss_D; // array of missings' positions

  //censoring
  real censL_val;
  int n_censL;                     // total number of obs at LB across D
  array[D] int n_censL_D;          // obs at LB per D
  array[D,max(n_censL_D)] int pos_censL_D; // array of obs at LBs' positions
  real censR_val;
  int n_censR;                      // total number of obs at LB across D
  array[D] int n_censR_D;           // obs at LB per D
  array[D,max(n_censR_D)] int pos_censR_D; // array of obs at LBs' positions

  // model adaptions based on user inputs:
  array[D_cen] int<lower=0, upper=1> innos_rand; // 1=person specific (random), 0=fixed
  int n_innos_fix;
  array[D_cen] int innos_pos;
  array[D_cen] int innos_fix_pos;


  // - dynamic model specification per D
  array[D] int<lower=0> N_pred;   // number of predictors per dimension
  array[D, max(N_pred)] int<lower=0> D_pred;   // matrix to index predictors to use per dimension
  array[D, max(N_pred)] int<lower=0> Lag_pred; // matrix to index lag of used predictors
  array[D] int Dpos1;  // index positions of dynamic effect parameters
  array[D] int Dpos2;
  array[D,max(N_pred)] int D_pred2;    // matrix to index predictors to use per dimension
  array[D,max(N_pred)] int Lag_pred2;  // matrix to index lag of used predictors

  // - time-invariant variables:
  // covariates as predictors of random effects
  int<lower=1> n_cov;           // number of covariates - minimum of 1 for intercepts
  int n_cov_bs;
  array[n_cov_bs, 2] int n_cov_mat;
  matrix[N, n_cov] W;  // predictors of individual parameters

  // outcome prediction
  int n_out;                 // number of outcome variables
  array[n_out,1] int n_out_bs;     // number of predictors per outcome
  int n_out_bs_max;          // number of predictors per outcome
  int n_out_bs_sum;          // number of predictors per outcome
  array[n_out,n_out_bs_max] int n_out_b_pos; // index positions
  int n_z;              // number of additional time-invariant as outcome predictors
  matrix[N, n_z] Z;     // observations of Z
  array[n_out] vector[N] out;        // outcome

  // group specific
  array[G] int N_G;  // number of clusters (persons) by group
  array[N] int g_id; // group index per cluster
  array[G, max(N_G)] int g_id_pos; // cluster (person) indexes by group

    // priors
  matrix[n_random,2] prior_gamma;
  matrix[n_random,2] prior_sd_R;
  real prior_LKJ;
  matrix[n_fixed,2] prior_b_fix;
//  matrix[n_fixed,2] prior_b_fix_diff;
  matrix[n_innos_fix,2] prior_sigma;
//  matrix[n_innos_fix,2] prior_sigma_diff;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;

  // covariances of innovations
  int n_inno_covs; // number of potential innovation covs to include
  int n_obs_cov;   // total number of residuals
  array[1,n_inno_covs] int inno_cov_pos;
  array[2] int<lower=-1,upper=1> inno_cov_load;

  array[D] int<lower=0,upper=1> is_wcen;   // parameter should be within centered = 1; should not = 0
  array[D] int<lower=0,upper=D> D_cen_pos; // pos of parameters that should be centered

  array[D] int<lower=0, upper=1> is_rdsem;
  array[D] int<lower=0> N_pred_rdsem;
  array[D, max(N_pred_rdsem)] int<lower=0> D_pred_rdsem;
  array[D] int Dpos1_rdsem;
  
  int<lower=1> grainsize;
}

transformed data{
  // creating pos and pos_cov for partial_sum
  array[N] int pos_start;
  array[N] int pos_end;
  array[N] int seq_N;
  array[N] int pos_start_cov;
  array[N] int pos_end_cov;
  int pos = 1;
  int pos_cov = 1;
  int obs_id_temp;
  int obs_id_temp_cov;

  for (n in 1:N){
    seq_N[n] = n;
    obs_id_temp = N_obs_id[n];
    obs_id_temp_cov = N_obs_id[n] - maxLag;

    pos_start[n] = pos;
    pos_start_cov[n] = pos_cov;
    pos_end[n] = pos + obs_id_temp -1;
    pos_end_cov[n] = pos_cov + obs_id_temp_cov - 1;

    pos = pos + obs_id_temp;
    pos_cov = pos_cov + obs_id_temp_cov;
  }

  // creating positioning for missings in partial sum
  array[N, D] int pos_start_miss = rep_array(0, N, D);  //starting position for each person and each dimension in y_impute
  array[N, D] int pos_end_miss = rep_array(0, N, D);    //ending postion
  array[N, D] int seq_N_miss = rep_array(0, N, D);      //number of missing values for each person and each dimension

  array[N, D] int pos_start_censL = rep_array(0, N, D);
  array[N, D] int pos_end_censL = rep_array(0, N, D);
  array[N, D] int seq_N_censL = rep_array(0, N, D);

  array[N, D] int pos_start_censR = rep_array(0, N, D);
  array[N, D] int pos_end_censR = rep_array(0, N, D);
  array[N, D] int seq_N_censR = rep_array(0, N, D);

  int impute_pos = 0;
  int impute_pos_censL = 0;
  int impute_pos_censR = 0;

  for (d in 1:D){
    for (n in 1:N){
      int counter = 0;
      int end_pos = 0;
      int counter_censL = 0;
      int end_pos_censL = 0;
      int counter_censR = 0;
      int end_pos_censR = 0;

      for (x in 1: n_miss_D[d]){
        if (pos_miss_D[d, x] >= pos_start[n] && pos_miss_D[d, x] <= pos_end[n]){
          if (counter == 0){
            pos_start_miss[n, d] = impute_pos + x;
          }
          counter = counter + 1;
          end_pos = x;
        }
      }
      for (x in 1:n_censL_D[d]){
        if (pos_censL_D[d, x] >= pos_start[n] && pos_censL_D[d,x] <= pos_end[n]){
          if(counter_censL == 0){
            pos_start_censL[n, d] = impute_pos_censL + x;
          }
          counter_censL = counter_censL + 1;
          end_pos_censL = x;
        }
      }
      for (x in 1:n_censR_D[d]){
        if (pos_censR_D[d, x] >= pos_start[n] && pos_censR_D[d,x] <= pos_end[n]){
          if(counter_censR == 0){
            pos_start_censR[n, d] = impute_pos_censR + x;
          }
          counter_censR = counter_censR + 1;
          end_pos_censR = x;
        }
      }
      if (counter > 0){
        pos_end_miss[n, d] = impute_pos + end_pos;
      }
      if (counter_censL > 0){
        pos_end_censL[n, d] = impute_pos_censL + end_pos_censL;
      }
      if (counter_censR > 0){
        pos_end_censR[n, d] = impute_pos_censR + end_pos_censR;
      }
      seq_N_miss[n, d] = counter;
      seq_N_censL[n, d] = counter_censL;
      seq_N_censR[n, d]  = counter_censR;
    }
    impute_pos = impute_pos + n_miss_D[d];
    impute_pos_censL = impute_pos_censL + n_censL_D[d];
    impute_pos_censR = impute_pos_censR + n_censR_D[d];
  }
}

parameters {
  array[N] vector[n_random] b_free;      // person-specific parameters
  array[G] vector[n_fixed] b_fix;        // fixed parameters
  array[G] vector<lower=0>[n_random] sd_R;        // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_mvn1] L;      // cholesky factor of random effects correlation matrix
  array[G] cholesky_factor_corr[n_mvn2] L2;     // cholesky factor of random effects correlation matrix
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
  array[G] vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  array[G] vector[n_out] alpha_out;               // outcome precition intercepts
  array[G] vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  array[G] vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  array[n_inno_covs] vector[n_obs_cov] eta_cov;
}

transformed parameters{
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  array[D_cen] vector[N] sd_noise;
  matrix[N, n_pars] b;
  array[n_inno_covs] vector[N] sd_inncov;

  // transformation of log-innovation variances if modeled as cluster-specific
  b = calculate_b(N, n_pars, n_random, is_random, b_free, n_fixed, G, g_id_pos,
  N_G, b_fix, is_fixed);

  bmu = calculate_bmu(G, N, n_random, gammas, n_cov, n_cov_bs, n_cov_mat, b_re_pred,
                     g_id_pos, W, N_G);


  for(d in 1:D_cen){
      if (innos_rand[d] == 0){
        for(g in 1:G){
          sd_noise[d, g_id_pos[g, 1:N_G[g]]] = rep_vector(sigma[g, innos_fix_pos[d]], N_G[g]);
        }
      }
      else{
        sd_noise[d] = sqrt(exp(b[,innos_pos[d]]));  // random effect transformed from log(var) to sd for each person
      }
  }
  // transform log innovation covarainces
  if(n_inno_covs > 0){
    for(i in 1:n_inno_covs){
      sd_inncov[i,1:N] = sqrt(exp(to_vector(b[,inno_cov_pos[1,i]])));
    }
  }
}


model {

  array[G] matrix[n_mvn1, n_mvn1] SIGMA;
  array[G] matrix[n_mvn2, n_mvn2] SIGMA2;

  for(g in 1:G){
    if(n_mvn1 > 0){
      SIGMA[g] = diag_pre_multiply(sd_R[g,pos_mvn1], L[g,]);
    }
    if(n_mvn2 > 0){
      SIGMA2[g] = diag_pre_multiply(sd_R[g,pos_mvn2], L2[g,]);
    }
  }

  // (Hyper-)Priors
  priors_lp(gammas, prior_gamma, sd_R, prior_sd_R, L, L2, prior_LKJ, sigma,
  n_innos_fix, prior_sigma, n_cov, b_re_pred, prior_b_re_pred, n_out, alpha_out, prior_alpha_out,
  b_out_pred, prior_b_out, sigma_out, prior_sigma_out, n_fixed, b_fix, prior_b_fix, n_mvn1, n_mvn2);

  
  target += reduce_sum(
    partial_sum_log_lik_kalman,
    seq_N,
    grainsize,
    N_obs_id, g_id, b_free, gammas, SIGMA, SIGMA2, D_cen, maxLag, D,
    is_wcen, y, pos_start, pos_end, pos_start_cov, pos_end_cov, b, D_cen_pos, N_pred,
    Lag_pred, D_pred, D_pred2, Lag_pred2, Dpos1, Dpos2, sd_noise, n_miss, n_miss_D,
    pos_miss_D, pos_start_miss, pos_end_miss, seq_N_miss, n_censL,
    n_censL_D, pos_censL_D, pos_start_censL, pos_end_censL,
    seq_N_censL, n_censR, n_censR_D, pos_censR_D, pos_start_censR,
    pos_end_censR, seq_N_censR, n_inno_covs, eta_cov, inno_cov_load, bmu, sd_R, sd_inncov, n_random,
    n_mvn1, n_mvn2, n_iid, pos_iid, pos_mvn1, pos_mvn2, is_rdsem, N_pred_rdsem, D_pred_rdsem, Dpos1_rdsem,
    censL_val, censR_val
  );



  // outcome prediction: get expectations of outcome values
  outcome_prediction_lp(n_out, G, n_random, n_z, N_G, g_id_pos, n_out_bs,
                        n_out_b_pos, b, Z, out, alpha_out, b_out_pred, sigma_out,
                        is_random);

}

generated quantities{
  array[G] matrix[n_mvn1,n_mvn1] bcorr;  // random coefficients correlation matrix
  array[G] matrix[n_mvn2,n_mvn2] bcorr2; // random coefficients correlation matrix
  for(g in 1:G){
    if(n_mvn1 > 0){
      bcorr[g] = multiply_lower_tri_self_transpose(L[g]);
    }
    if(n_mvn2 > 0){
      bcorr2[g] = multiply_lower_tri_self_transpose(L2[g]);
    }
    }
}

