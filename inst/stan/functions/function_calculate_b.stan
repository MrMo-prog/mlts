matrix calculate_b(int N,
                  int n_pars,
                  int n_random,
                  array[] int is_random,
                  array[] vector b_free,
                  int n_fixed,
                  int G,
                  array[,] int g_id_pos,
                  array[] int N_G,
                  array[] vector b_fix,
                  array[,] int is_fixed){

  matrix[N, n_pars] b;
  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      for(g in 1:G){
        b[g_id_pos[g,1:N_G[g]],is_fixed[1,i]] = rep_vector(b_fix[g,i],N_G[g]);
      }
    }
  }
  return b;
}
