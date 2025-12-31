void outcome_prediction_lp(int n_out, int G, int n_random, int n_z,
                           array[] int N_G, array[,] int g_id_pos,
                           array[,] int n_out_bs, array[,] int n_out_b_pos,
                           matrix b, matrix Z, array[] vector out,
                           array[] vector alpha_out, array[] vector b_out_pred,
                           array[] vector sigma_out, array[] int is_random){

  if(n_out > 0){
    for(g in 1:G){
      int k = 1;
      matrix[N_G[g], n_random + n_z] b_z = append_col(b[g_id_pos[g, 1:N_G[g]], is_random], Z[g_id_pos[g, 1:N_G[g]], ]);
      for(i in 1:n_out){
        int n_bs = n_out_bs[i,1];
        target += normal_lpdf(out[i, g_id_pos[g, 1:N_G[g]]] | alpha_out[g,i] + b_z[, n_out_b_pos[i, 1:n_bs]]
                              * segment(b_out_pred[g,], k, n_bs), sigma_out[g,i]);
        k = k + n_bs;
      }
    }
  }
}
