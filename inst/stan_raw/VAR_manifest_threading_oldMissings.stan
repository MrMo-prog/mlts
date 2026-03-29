functions{
  #include "../stan/functions/threading/partial_sum_without_missings.stan"
  #include "../stan/functions/function_calculate_b.stan"
  #include "../stan/functions/function_hyper_priors.stan"
  #include "../stan/functions/function_calculate_bmu.stan"
  #include "../stan/functions/function_outcome_prediction.stan"
  #include "../stan/functions/function_missings_and_censoring.stan"
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
}

parameters {
  array[N] vector[n_random] b_free;      // person-specific parameters
  array[G] vector[n_fixed] b_fix;        // fixed parameters
  array[G] vector<lower=0>[n_random] sd_R;        // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_random] L;      // cholesky factor of random effects correlation matrix
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
  array[G] vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  array[G] vector[n_out] alpha_out;               // outcome precition intercepts
  array[G] vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  array[G] vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  vector[n_miss] y_impute;                        // vector to store imputed values
  array[n_inno_covs] vector[n_obs_cov] eta_cov;
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<lower=censR_val>[n_censR] y_impute_censR;
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
  array[D] vector[N_obs] y_merge = y;
  array[G] matrix[n_random, n_random] SIGMA;

  for(g in 1:G){
    SIGMA[g] = diag_pre_multiply(sd_R[g], L[g]); // covariance matrix of parameters by group
  }

  // missings and censoring
  if (n_miss > 0){
    y_merge = missings_and_censoring(y_merge, n_miss_D, pos_miss_D, y_impute);
  }
  if (n_censL > 0){
    y_merge = missings_and_censoring(y_merge, n_censL_D, pos_censL_D, y_impute_censL);
  }
  if (n_censR > 0){
    y_merge = missings_and_censoring(y_merge, n_censR_D, pos_censR_D, y_impute_censR);
  }

  // (Hyper-)Priors
  priors_lp(gammas, prior_gamma, sd_R, prior_sd_R, L, prior_LKJ, sigma,
  n_innos_fix, prior_sigma, n_cov, b_re_pred, prior_b_re_pred, n_out, alpha_out, prior_alpha_out,
  b_out_pred, prior_b_out, sigma_out, prior_sigma_out, n_fixed, b_fix, prior_b_fix);

  target += reduce_sum(
    partial_sum_log_lik,
    seq_N,
    grainsize,
    N_obs_id, g_id, b_free, gammas, SIGMA, D_cen, maxLag, D,
    is_wcen, y_merge, pos_start, pos_end, pos_start_cov, pos_end_cov, b, D_cen_pos, N_pred,
    Lag_pred, D_pred, D_pred2, Lag_pred2, Dpos1, Dpos2, sd_noise, n_inno_covs,
    eta_cov, inno_cov_load, bmu, sd_R, sd_inncov, n_random
  );


  // outcome prediction: get expectations of outcome values
  outcome_prediction_lp(n_out, G, n_random, n_z, N_G, g_id_pos, n_out_bs,
                        n_out_b_pos, b, Z, out, alpha_out, b_out_pred, sigma_out,
                        is_random);

}

generated quantities{
  array[G] matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
    for(g in 1:G){
        bcorr[g] = multiply_lower_tri_self_transpose(L[g]);
      }
}

