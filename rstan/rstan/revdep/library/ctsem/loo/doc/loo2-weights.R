params <-
list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

## ---- setup, message=FALSE-----------------------------------------------
library(rstanarm)
library(loo)

## ------------------------------------------------------------------------
data(milk)
d <- milk[complete.cases(milk),]
d$neocortex <- d$neocortex.perc /100
str(d)

## ---- results="hide"-----------------------------------------------------
fit1 <- stan_glm(kcal.per.g ~ 1, data = d, seed = 2030)
fit2 <- update(fit1, formula = kcal.per.g ~ neocortex)
fit3 <- update(fit1, formula = kcal.per.g ~ log(mass))
fit4 <- update(fit1, formula = kcal.per.g ~ neocortex + log(mass))

## ------------------------------------------------------------------------
waic1 <- waic(fit1)
waic2 <- waic(fit2)
waic3 <- waic(fit3)
waic4 <- waic(fit4)
waics <- c(
  waic1$estimates["elpd_waic", 1],
  waic2$estimates["elpd_waic", 1],
  waic3$estimates["elpd_waic", 1],
  waic4$estimates["elpd_waic", 1]
)

## ------------------------------------------------------------------------
# note: the loo function accepts a 'cores' argument that we recommend specifying
# when working with bigger datasets

loo1 <- loo(fit1)
loo2 <- loo(fit2)
loo3 <- loo(fit3)
loo4 <- loo(fit4)
lpd_point <- cbind(
  loo1$pointwise[,"elpd_loo"], 
  loo2$pointwise[,"elpd_loo"],
  loo3$pointwise[,"elpd_loo"], 
  loo4$pointwise[,"elpd_loo"]
)

## ------------------------------------------------------------------------
print(loo3)
print(loo4)

## ------------------------------------------------------------------------
waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2)

## ------------------------------------------------------------------------
waic_wts_demo <- 
  exp(waics[c(1,1,1,1,1,1,1,1,1,1,2,3,4)]) /
  sum(exp(waics[c(1,1,1,1,1,1,1,1,1,1,2,3,4)]))
round(waic_wts_demo, 3)

## ------------------------------------------------------------------------
stacking_weights(lpd_point[,c(1,1,1,1,1,1,1,1,1,1,2,3,4)])

## ------------------------------------------------------------------------
data(Kline)
d <- Kline
d$log_pop <- log(d$population)
d$contact_high <- ifelse(d$contact=="high", 1, 0)
str(d)

## ---- results="hide"-----------------------------------------------------
fit10 <-
  stan_glm(
    total_tools ~ log_pop + contact_high + log_pop * contact_high,
    family = poisson(link = "log"),
    data = d,
    prior = normal(0, 1, autoscale = FALSE),
    prior_intercept = normal(0, 100, autoscale = FALSE),
    seed = 2030
  )

## ------------------------------------------------------------------------
loo10 <- loo(fit10)
print(loo10)

## ------------------------------------------------------------------------
loo10 <- loo(fit10, k_threshold=0.7)
print(loo10)

## ------------------------------------------------------------------------
waic10 <- waic(fit10)
print(waic10)

## ---- results="hide"-----------------------------------------------------
fit11 <- update(fit10, formula = total_tools ~ log_pop + contact_high)
fit12 <- update(fit10, formula = total_tools ~ log_pop)

## ------------------------------------------------------------------------
(loo11 <- loo(fit11))
(loo12 <- loo(fit12))

## ------------------------------------------------------------------------
loo11 <- loo(fit11, k_threshold=0.7)
loo12 <- loo(fit12, k_threshold=0.7)
lpd_point <- cbind(
  loo10$pointwise[, "elpd_loo"], 
  loo11$pointwise[, "elpd_loo"], 
  loo12$pointwise[, "elpd_loo"]
)

## ------------------------------------------------------------------------
waic11 <- waic(fit11)
waic12 <- waic(fit12)
waics <- c(
  waic10$estimates["elpd_waic", 1], 
  waic11$estimates["elpd_waic", 1], 
  waic12$estimates["elpd_waic", 1]
)

## ------------------------------------------------------------------------
waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2)

## ------------------------------------------------------------------------
# using list of loo objects
loo_list <- list(loo10, loo11, loo12)
loo_model_weights(loo_list)
loo_model_weights(loo_list, method = "pseudobma")
loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)

