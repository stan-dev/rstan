## this example is mainly for Windows since function mclapply
## in package parallel does not work on Windows 

library(parallel)
library(rstan)

foo_data <- list(N = 3)
rng_seed <- 234
foo <- stan(model_code = 'data { int N; }  parameters { vector[N] y; } model { y ~ normal(0,1); }',
            data = foo_data, chains = 0)

CL = makeCluster(2)
clusterExport(cl = CL, c("foo_data", "foo", "rng_seed")) 
sflist1 <- parLapply(CL, 1:4, fun = function(cid) {  
  require(rstan)
  stan(fit = foo, data = foo_data, chains = 1, iter = 2000, seed = rng_seed, chain_id = cid)
})

fit <- sflist2stanfit(sflist1)
print(fit)
stopCluster(CL)
