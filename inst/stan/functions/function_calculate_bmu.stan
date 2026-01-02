matrix calculate_bmu(int G,
                    int N,
                    int n_random,
                    array[] row_vector gammas,
                    int n_cov,
                    int n_cov_bs,
                    array[,] int n_cov_mat,
                    array[] vector b_re_pred,
                    array[,] int g_id_pos,
                    matrix W,
                    array[] int N_G){

  matrix[N, n_random] bmu;
  array[G] matrix[n_cov, n_random] b_re_pred_mat;

   // REs regressed on covariates
  for(g in 1:G){
    b_re_pred_mat[g] = rep_matrix(0, n_cov, n_random);
    b_re_pred_mat[g,1,] = gammas[g,];
    if(n_cov>1){
      for(i in 1:n_cov_bs){
      b_re_pred_mat[g,n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[g,i];
      }
    }
    // calculate population means (intercepts) of person-specific parameters
    bmu[g_id_pos[g,1:N_G[g]],] = W[g_id_pos[g,1:N_G[g]],] * b_re_pred_mat[g];
  }
  return bmu;
}
