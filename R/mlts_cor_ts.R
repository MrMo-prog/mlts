mlts_cor_ts <- function(data, ts, id, type = "between", as_mat = TRUE, within_fn = NULL){
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
            c_test <- tryCatch({
              stats::cor.test(person_means[, ts[item]],
                              person_means[, ts[item2]], method = "pearson")
            }, error = function(e) {
              list(conf.int = c(NA, NA))
            })


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
  else if (type == "within"){
    u_ids <- unique(data[[id]])
    matrices <- list()
    for ( i in u_ids){
      person_data <- data[data[[id]] == i, ts]
      matrices[[as.character(i)]] <- cor(person_data, method = "pearson", use = "pairwise.complete.obs")
    }
    if(as_mat){
      if (!is.null(within_fn)){
        temp_cors <- numeric(length(matrices))
        final_mat <- matrix(nrow = n_items, ncol = n_items)
        colnames(final_mat) <- ts
        rownames(final_mat) <- ts

        for (item in 1:(n_items-1)){
          for (item2 in (item+1):n_items){
            for ( mat in seq_along(matrices)){
              temp <- matrices[[mat]]
              temp_cors[mat] <- temp[item, item2]
            }
            fn_value <- within_fn(temp_cors, na.rm = TRUE)
            final_mat[item, item2] <- fn_value
            final_mat[item2, item] <- fn_value
          }
        }
        diag(final_mat) <- 1
        return(final_mat)
      }
      return(matrices)
    }
    else {
      results <- list()
      list_counter = 1
      for (i in u_ids){
        person_data <- data[data[[id]] == i, ts]
        for (item in 1:(n_items-1)){
          for (item2 in (item+1):n_items){
            vec1 <- person_data[, item]
            vec2 <- person_data[, item2]

            c_test <- tryCatch({
                stats::cor.test(vec1,vec2, method = "pearson")
              }, error = function(e){
                  list(estimate = NA, conf.int = c(NA, NA))
                })

            temp_df <- data.frame(
              ID = i,
              Var1 = ts[item],
              Var2 = ts[item2],
              r = c_test$estimate,
              CI.lb = c_test$conf.int[1],
              CI.ub = c_test$conf.int[2]
            )
            results[[list_counter]] <- temp_df
            list_counter <- list_counter +1
          }
        }
      }
      return(do.call(rbind, results))
    }
  }

}
