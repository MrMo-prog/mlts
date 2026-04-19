mlts_rel <- function(fit_data, method=1){

  if (method==1){
    split_data <- split(fit_data$person.pars.summary, fit_data$person.pars.summary$Param)
    reliability_results <- sapply(split_data, function(x) {
      var_EAP <- var(x$mean, na.rm = TRUE)
      mean_err_var <- mean(x$sd^2, na.rm = TRUE)
      rel <- var_EAP / (var_EAP + mean_err_var)
      return(rel)
    })
    return(reliability_results)
  }
  else {
    # new method
    draws_rstan <- rstan::extract(fit_data$stanfit)
    b_free_array <- draws_rstan$b_free

    # 1. preparation of a list of all N*K matrices
    name_vec <- fit_data$model$Param[fit_data$model$Level == "Within" & fit_data$model$isRandom==1]
    matrix_list <- asplit(b_free_array, MARGIN = 3)
    matrix_list <- lapply(matrix_list, t) # list of all N*K matrices
    names(matrix_list) <- name_vec

    # 2. randomly sort columns and split matrices into M and W
    matrix_list <- lapply(matrix_list, function(mat){
      mat[, sample(ncol(mat))]
    })
    split_matrix_list <- lapply(matrix_list, function (mat){
      mid <- ncol(mat) / 2
      list(
        M = mat[, 1:mid],
        W = mat[, (mid+1):ncol(mat)]
      )
    })

    # 3. generate vector of pearson correlations between M and W
    cor_list <- lapply(split_matrix_list, function (mat){
      cor_vec <- sapply(1:ncol(mat[[1]]), function(i){
        cor_val <- stats::cor(mat$M[, i], mat$W[, i], method = "pearson")
        return(cor_val)
      })
    })
    # 4. calculating reliability etc. from p
    reliability_results <- lapply(cor_list, function(i){
      data.frame (
        Mean_Rel = mean(i),
        SD_Rel = sd(i)
      )
    })
    return(reliability_results)
  }

}


