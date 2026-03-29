cor_ts <- function(data, ts, id, type = "between", as_mat = TRUE){
  if (length(ts) == 1){
    stop("ts must contain at least two items to calculate a correlation", call. = FALSE)
  }
  n_items <- length(ts)

  # between level correlations
  if(type == "between"){
    person_means <- aggregate(x = data[, ts ], by = list(data[[id]]), FUN = mean, na.rm = TRUE)
    cor_mat <- cor(person_means[,-1], method = "pearson", use = "pairwise.complete.obs")
    if (as_mat){
      return(cor_mat)
    } else {
      results <- list()
      list_counter <- 1
      for (item in 1:(n_items-1)){
        for (item2 in (item+1): n_items){
            c_test <- stats::cor.test(person_means[, ts[item]],
                                      person_means[, ts[item2]], method = "pearson")
            temp_df <- data.frame(
              Var1 = ts[item],
              Var2 = ts[item2],
              r = cor_mat[item, item2],
              CI.lb = c_test$conf.int[1],
              CI.ub = c_test$conf.int[2]
            )
            results[[list_counter]] <- temp_df
            list_counter <- list_counter + 1
        }
      }
      return(do.call(rbind, results))
    }
  }

}
