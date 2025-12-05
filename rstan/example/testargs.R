library(rstan)
sm <- stan_model(model_code = 'data { int N; } parameters { real y[N]; } model {y ~ normal(0, 1); }')

a <- sampling(sm, data = list(N = 3))

b <- sampling(sm, controls = list(metric = "d"))

c <- sampling(sm, control = list(metric = "e"))

c2 <- sampling(sm, data = list(N = 3), control = list(metric = "e"))

c3 <- sampling(sm, data = list(N = 3), control = list(unmetric = "e"))

d <- stan(fit = a, unkown_arg1 = 1, unkown_arg2 = 2, data = list(N = 4))
