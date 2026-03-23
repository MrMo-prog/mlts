library(mlts)
m1 <- mlts_model(q = 1)

simData <- mlts_sim(TP=400, N=200, model=m1, seed=113, default=TRUE)
fitted <- mlts_fit(model=m1, data=simData, id = "ID", ts = "Y1", iter=1000)
fitted_threading <- mlts_fit(model=m1, data=simData, id = "ID", ts = "Y1", iter=1000, threads_per_chain = 2, grainsize = 2)

apply(rstan::get_elapsed_time(fitted$stanfit),1,sum)
apply(rstan::get_elapsed_time(fitted_threading$stanfit),1,sum)
1 - ( max(apply(rstan::get_elapsed_time(fitted$stanfit),1,sum)) /
        max(apply(rstan::get_elapsed_time(fitted_threading$stanfit),1,sum)) )

apply(round(rstan::summary(fitted_threading$stanfit)$summary,1) == round(rstan::summary(fitted$stanfit)$summary,1),2,
      function(x){sum(x,na.rm = T)/sum(!is.na(x))})
