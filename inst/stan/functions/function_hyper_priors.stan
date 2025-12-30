void priors_lp (array[] row_vector gammas,
                matrix prior_gamma,
                array[] vector sd_R,
                matrix prior_sd_R,
                array[] matrix L,
                real prior_LKJ,

                array[] vector sigma,
                int n_innos_fix,
                matrix prior_sigma,

                int n_cov,
                array[] vector b_re_pred,
                matrix prior_b_re_pred,

                int n_out,
                array[] vector alpha_out,
                matrix prior_alpha_out,
                array[] vector b_out_pred,
                matrix prior_b_out,
                array[] vector sigma_out,
                matrix prior_sigma_out,

                int n_fixed,
                array[] vector b_fix,
                matrix prior_b_fix){

  int G = size(gammas);

  for(g in 1:G){
    target += normal_lpdf(gammas[g] | prior_gamma[,1],prior_gamma[,2]);
    target += cauchy_lpdf(sd_R[g] | prior_sd_R[,1], prior_sd_R[,2]);
    target += lkj_corr_cholesky_lpdf(L[g] | prior_LKJ);

    if(n_innos_fix>0){
      target += cauchy_lpdf(sigma[g] | prior_sigma[,1], prior_sigma[,2]);
    }
    if(n_cov > 1){
      target += normal_lpdf(b_re_pred[g] | prior_b_re_pred[,1], prior_b_re_pred[,2]);
    }
    if(n_out > 0){
      target += normal_lpdf(alpha_out[g] | prior_alpha_out[,1], prior_alpha_out[,2]);
      target += normal_lpdf(b_out_pred[g] | prior_b_out[,1], prior_b_out[,2]);
      target += cauchy_lpdf(sigma_out[g] | prior_sigma_out[,1], prior_sigma_out[,2]);
    }
    if(n_fixed > 0){
      target += normal_lpdf(b_fix[g] | prior_b_fix[,1],prior_b_fix[,2]);
      }
    }
}

// zusätzlicher Parameter L_inno für covsfix Modelle (gleiche Funktion)

void priors_lp (array[] row_vector gammas,
                matrix prior_gamma,
                array[] vector sd_R,
                matrix prior_sd_R,
                array[] matrix L,
                real prior_LKJ,

                array[] vector sigma,
                int n_innos_fix,
                matrix prior_sigma,

                int n_cov,
                array[] vector b_re_pred,
                matrix prior_b_re_pred,

                int n_out,
                array[] vector alpha_out,
                matrix prior_alpha_out,
                array[] vector b_out_pred,
                matrix prior_b_out,
                array[] vector sigma_out,
                matrix prior_sigma_out,

                int n_fixed,
                array[] vector b_fix,
                matrix prior_b_fix,

                array[] matrix L_inno //nur angeben bei covsfix Modellen
                ){
  // Basisfunktion
  priors_lp(gammas, prior_gamma, sd_R, prior_sd_R, L, prior_LKJ,
            sigma, n_innos_fix, prior_sigma, n_cov, b_re_pred, prior_b_re_pred,
            n_out, alpha_out, prior_alpha_out, b_out_pred, prior_b_out, sigma_out, prior_sigma_out,
            n_fixed, b_fix, prior_b_fix);

  // neue Variable für covsfix Modelle
  int G = size(gammas);
  for(g in 1:G){
    target += lkj_corr_cholesky_lpdf(L_inno[g] | prior_LKJ);
    }
}
