pkgname <- "broom"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('broom')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("argument_glossary")
### * argument_glossary

flush(stderr()); flush(stdout())

### Name: argument_glossary
### Title: Allowed argument names in tidiers
### Aliases: argument_glossary
### Keywords: datasets

### ** Examples

argument_glossary



cleanEx()
nameEx("augment.decomposed.ts")
### * augment.decomposed.ts

flush(stderr()); flush(stdout())

### Name: augment.decomposed.ts
### Title: Augment data with information from a(n) decomposed.ts object
### Aliases: augment.decomposed.ts decompose_tidiers

### ** Examples


# Time series of temperatures in Nottingham, 1920-1939:
nottem

# Perform seasonal decomposition on the data with both decompose
# and stl:
d1 <- stats::decompose(nottem)
d2 <- stats::stl(nottem, s.window = "periodic", robust = TRUE)

# Compare the original series to its decompositions.

cbind(broom::tidy(nottem), broom::augment(d1),
      broom::augment(d2))

# Visually compare seasonal decompositions in tidy data frames.

library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

decomps <- tibble(
    # Turn the ts objects into data frames.
    series = list(as.data.frame(nottem), as.data.frame(nottem)),
    # Add the models in, one for each row.
    decomp = c("decompose", "stl"),
    model = list(d1, d2)
) %>%
    rowwise() %>%
    # Pull out the fitted data using broom::augment.
    mutate(augment = list(broom::augment(model))) %>%
    ungroup() %>%
    # Unnest the data frames into a tidy arrangement of
    # the series next to its seasonal decomposition, grouped
    # by the method (stl or decompose).
    group_by(decomp) %>%
    unnest(series, augment) %>%
    mutate(index = 1:n()) %>%
    ungroup() %>%
    select(decomp, index, x, adjusted = .seasadj)

ggplot(decomps) +
    geom_line(aes(x = index, y = x), colour = "black") +
    geom_line(aes(x = index, y = adjusted, colour = decomp,
                  group = decomp))




cleanEx()
nameEx("augment.ivreg")
### * augment.ivreg

flush(stderr()); flush(stdout())

### Name: augment.ivreg
### Title: Augment data with information from a(n) ivreg object
### Aliases: augment.ivreg

### ** Examples


library(AER)

data("CigarettesSW", package = "AER")
ivr <- ivreg(
  log(packs) ~ income | population,
  data = CigarettesSW,
  subset = year == "1995"
)

summary(ivr)

tidy(ivr)
tidy(ivr, conf.int = TRUE)
tidy(ivr, conf.int = TRUE, exponentiate = TRUE)

augment(ivr)

glance(ivr)




cleanEx()
nameEx("augment.loess")
### * augment.loess

flush(stderr()); flush(stdout())

### Name: augment.loess
### Title: Tidy a(n) loess object
### Aliases: augment.loess loess_tidiers

### ** Examples


lo <- loess(mpg ~ wt, mtcars)
augment(lo)

# with all columns of original data
augment(lo, mtcars)

# with a new dataset
augment(lo, newdata = head(mtcars))




cleanEx()
nameEx("augment.smooth.spline")
### * augment.smooth.spline

flush(stderr()); flush(stdout())

### Name: augment.smooth.spline
### Title: Tidy a(n) smooth.spline object
### Aliases: augment.smooth.spline smooth.spline_tidiers

### ** Examples


spl <- smooth.spline(mtcars$wt, mtcars$mpg, df = 4)
augment(spl, mtcars)
augment(spl)  # calls original columns x and y

library(ggplot2)
ggplot(augment(spl, mtcars), aes(wt, mpg)) +
    geom_point() + geom_line(aes(y = .fitted))




cleanEx()
nameEx("bootstrap")
### * bootstrap

flush(stderr()); flush(stdout())

### Name: bootstrap
### Title: Set up bootstrap replicates of a dplyr operation
### Aliases: bootstrap

### ** Examples


## Not run: 
##D library(dplyr)
##D mtcars %>% bootstrap(10) %>% do(tidy(lm(mpg ~ wt, .)))
## End(Not run)




cleanEx()
nameEx("brms_tidiers")
### * brms_tidiers

flush(stderr()); flush(stdout())

### Name: brms_tidiers
### Title: Tidying methods for a brms model
### Aliases: brms_tidiers tidy.brmsfit

### ** Examples

## Not run: 
##D  library(brms)
##D  fit <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
##D             iter = 500, chains = 2)
##D  tidy(fit)
##D  tidy(fit, parameters = "^sd_", intervals = FALSE)
##D  tidy(fit, par_type = "non-varying")
##D  tidy(fit, par_type = "varying")
##D  tidy(fit, par_type = "hierarchical", robust = TRUE)
## End(Not run)




cleanEx()
nameEx("column_glossary")
### * column_glossary

flush(stderr()); flush(stdout())

### Name: column_glossary
### Title: Allowed column names in tidied tibbles
### Aliases: column_glossary
### Keywords: datasets

### ** Examples

column_glossary



cleanEx()
nameEx("data.frame_tidiers")
### * data.frame_tidiers

flush(stderr()); flush(stdout())

### Name: data.frame_tidiers
### Title: Tidiers for data.frame objects
### Aliases: data.frame_tidiers tidy.data.frame glance.data.frame

### ** Examples


## Not run: 
##D td <- tidy(mtcars)
##D td
##D 
##D glance(mtcars)
##D 
##D library(ggplot2)
##D # compare mean and standard deviation
##D ggplot(td, aes(mean, sd)) + geom_point() +
##D      geom_text(aes(label = column), hjust = 1, vjust = 1) +
##D      scale_x_log10() + scale_y_log10() + geom_abline()
## End(Not run)




cleanEx()
nameEx("durbinWatsonTest_tidiers")
### * durbinWatsonTest_tidiers

flush(stderr()); flush(stdout())

### Name: durbinWatsonTest_tidiers
### Title: Tidy/glance a(n) durbinWatsonTest object
### Aliases: durbinWatsonTest_tidiers tidy.durbinWatsonTest
###   glance.durbinWatsonTest

### ** Examples


dw <- car::durbinWatsonTest(lm(mpg ~ wt, data = mtcars))
tidy(dw)
glance(dw)  # same output for all durbinWatsonTests




cleanEx()
nameEx("emmeans_tidiers")
### * emmeans_tidiers

flush(stderr()); flush(stdout())

### Name: emmeans_tidiers
### Title: Tidy estimated marginal means (least-squares means) objects from
###   the emmeans and lsmeans packages
### Aliases: emmeans_tidiers tidy.lsmobj tidy.ref.grid tidy.emmGrid

### ** Examples


if (require("emmeans", quietly = TRUE)) {
  # linear model for sales of oranges per day
  oranges_lm1 <- lm(sales1 ~ price1 + price2 + day + store, data = oranges)

  # reference grid; see vignette("basics", package = "emmeans")
  oranges_rg1 <- ref_grid(oranges_lm1)
  td <- tidy(oranges_rg1)
  td

  # marginal averages
  marginal <- emmeans(oranges_rg1, "day")
  tidy(marginal)

  # contrasts
  tidy(contrast(marginal))
  tidy(contrast(marginal, method = "pairwise"))

  # plot confidence intervals
  library(ggplot2)
  ggplot(tidy(marginal), aes(day, estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high))

  # by multiple prices
  by_price <- emmeans(oranges_lm1, "day", by = "price2",
                      at = list(price1 = 50, price2 = c(40, 60, 80),
                      day = c("2", "3", "4")) )
  by_price
  tidy(by_price)

  ggplot(tidy(by_price), aes(price2, estimate, color = day)) +
    geom_line() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high))
}




cleanEx()
nameEx("geeglm_tidiers")
### * geeglm_tidiers

flush(stderr()); flush(stdout())

### Name: tidy.geeglm
### Title: Tidy a(n) geeglm object
### Aliases: tidy.geeglm geeglm_tidiers geepack_tidiers

### ** Examples


if (requireNamespace("geepack", quietly = TRUE)) {
  library(geepack)
  data(state)
    
  ds <- data.frame(state.region, state.x77)

  geefit <- geeglm(Income ~ Frost + Murder, id = state.region,
                   data = ds, family = gaussian,
                   corstr = "exchangeable")

  tidy(geefit)
  tidy(geefit, quick = TRUE)
  tidy(geefit, conf.int = TRUE)
}




cleanEx()
nameEx("glance.betareg")
### * glance.betareg

flush(stderr()); flush(stdout())

### Name: glance.betareg
### Title: Glance at a(n) betareg object
### Aliases: glance.betareg

### ** Examples


library(betareg)

data("GasolineYield", package = "betareg")

mod <- betareg(yield ~ batch + temp, data = GasolineYield)

mod
tidy(mod)
tidy(mod, conf.int = TRUE)
tidy(mod, conf.int = TRUE, conf.level = .99)

augment(mod)

glance(mod)




cleanEx()
nameEx("glance.binDesign")
### * glance.binDesign

flush(stderr()); flush(stdout())

### Name: glance.binDesign
### Title: Glance at a(n) binDesign object
### Aliases: glance.binDesign

### ** Examples


if (require("binGroup", quietly = TRUE)) {
    des <- binDesign(nmax = 300, delta = 0.06,
                     p.hyp = 0.1, power = .8)

    glance(des)
    tidy(des)

    # the ggplot2 equivalent of plot(des)
    library(ggplot2)
    ggplot(tidy(des), aes(n, power)) +
        geom_line()
}




cleanEx()
nameEx("glance.glm")
### * glance.glm

flush(stderr()); flush(stdout())

### Name: glance.glm
### Title: Glance at a(n) glm object
### Aliases: glance.glm

### ** Examples


g <- glm(am ~ mpg, mtcars, family = "binomial")
glance(g)




cleanEx()
nameEx("glance.ivreg")
### * glance.ivreg

flush(stderr()); flush(stdout())

### Name: glance.ivreg
### Title: Glance at a(n) ivreg object
### Aliases: glance.ivreg

### ** Examples


library(AER)

data("CigarettesSW", package = "AER")
ivr <- ivreg(
  log(packs) ~ income | population,
  data = CigarettesSW,
  subset = year == "1995"
)

summary(ivr)

tidy(ivr)
tidy(ivr, conf.int = TRUE)
tidy(ivr, conf.int = TRUE, exponentiate = TRUE)

augment(ivr)

glance(ivr)




cleanEx()
nameEx("glance.lavaan")
### * glance.lavaan

flush(stderr()); flush(stdout())

### Name: glance.lavaan
### Title: Glance at a(n) lavaan object
### Aliases: glance.lavaan

### ** Examples


if (require("lavaan", quietly = TRUE)) {

 library(lavaan)

 cfa.fit <- cfa(
   'F =~ x1 + x2 + x3 + x4 + x5',
   data = HolzingerSwineford1939, group = "school"
 )
 glance(cfa.fit)

}




cleanEx()
nameEx("glance.rlm")
### * glance.rlm

flush(stderr()); flush(stdout())

### Name: glance.rlm
### Title: Glance at a(n) rlm object
### Aliases: glance.rlm rlm_tidiers

### ** Examples


library(MASS)

r <- rlm(stack.loss ~ ., stackloss)
tidy(r)
augment(r)
glance(r)




cleanEx()
nameEx("lme4_tidiers")
### * lme4_tidiers

flush(stderr()); flush(stdout())

### Name: lme4_tidiers
### Title: Tidying methods for mixed effects models
### Aliases: lme4_tidiers tidy.merMod augment.merMod glance.merMod

### ** Examples


## Not run: 
##D if (require("lme4")) {
##D     # example regressions are from lme4 documentation
##D     lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##D     tidy(lmm1)
##D     tidy(lmm1, effects = "fixed")
##D     tidy(lmm1, effects = "fixed", conf.int=TRUE)
##D     tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="profile")
##D     tidy(lmm1, effects = "ran_modes", conf.int=TRUE)
##D     head(augment(lmm1, sleepstudy))
##D     glance(lmm1)
##D 
##D     glmm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##D                   data = cbpp, family = binomial)
##D     tidy(glmm1)
##D     tidy(glmm1, effects = "fixed")
##D     head(augment(glmm1, cbpp))
##D     glance(glmm1)
##D 
##D     startvec <- c(Asym = 200, xmid = 725, scal = 350)
##D     nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
##D                   Orange, start = startvec)
##D     tidy(nm1)
##D     tidy(nm1, effects = "fixed")
##D     head(augment(nm1, Orange))
##D     glance(nm1)
##D }
## End(Not run)



cleanEx()
nameEx("matrix_tidiers")
### * matrix_tidiers

flush(stderr()); flush(stdout())

### Name: matrix_tidiers
### Title: Tidiers for matrix objects
### Aliases: matrix_tidiers tidy.matrix glance.matrix

### ** Examples


## Not run: 
##D mat <- as.matrix(mtcars)
##D tidy(mat)
##D glance(mat)
## End(Not run)




cleanEx()
nameEx("mcmc_tidiers")
### * mcmc_tidiers

flush(stderr()); flush(stdout())

### Name: mcmc_tidiers
### Title: Tidying methods for MCMC (Stan, JAGS, etc.) fits
### Aliases: mcmc_tidiers tidyMCMC tidy.rjags tidy.stanfit

### ** Examples


## Not run: 
##D 
##D # Using example from "RStan Getting Started"
##D # https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
##D 
##D model_file <- system.file("extdata", "8schools.stan", package = "broom")
##D 
##D schools_dat <- list(J = 8,
##D                     y = c(28,  8, -3,  7, -1,  1, 18, 12),
##D                     sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
##D 
##D if (requireNamespace("rstan", quietly = TRUE)) {
##D   set.seed(2015)
##D   rstan_example <- stan(file = model_file, data = schools_dat,
##D                         iter = 100, chains = 2)
##D }
##D 
## End(Not run)

if (requireNamespace("rstan", quietly = TRUE)) {
  # the object from the above code was saved as rstan_example.rda
  infile <- system.file("extdata", "rstan_example.rda", package = "broom")
  load(infile)

  tidy(rstan_example)
  tidy(rstan_example, conf.int = TRUE, pars = "theta")

  td_mean <- tidy(rstan_example, conf.int = TRUE)
  td_median <- tidy(rstan_example, conf.int = TRUE, estimate.method = "median")

  library(dplyr)
  library(ggplot2)
  tds <- rbind(mutate(td_mean, method = "mean"),
               mutate(td_median, method = "median"))

  ggplot(tds, aes(estimate, term)) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
    geom_point(aes(color = method))
}





cleanEx()
nameEx("mgcv_tidy_gam")
### * mgcv_tidy_gam

flush(stderr()); flush(stdout())

### Name: tidy.gam
### Title: Tidy a(n) gam object
### Aliases: tidy.gam mgcv_tidiers gam_tidiers

### ** Examples


g <- mgcv::gam(mpg ~ s(hp) + am + qsec, data = mtcars)
  
tidy(g)
tidy(g, parametric = TRUE)
glance(g)





cleanEx()
nameEx("nlme_tidiers")
### * nlme_tidiers

flush(stderr()); flush(stdout())

### Name: nlme_tidiers
### Title: Tidying methods for mixed effects models
### Aliases: nlme_tidiers tidy.lme augment.lme glance.lme

### ** Examples


## Not run: 
##D if (require("nlme") & require("lme4")) {
##D     # example regressions are from lme4 documentation, but used for nlme
##D     lmm1 <- lme(Reaction ~ Days, random=~ Days|Subject, sleepstudy)
##D     tidy(lmm1)
##D     tidy(lmm1, effects = "fixed")
##D     head(augment(lmm1, sleepstudy))
##D     glance(lmm1)
##D 
##D 
##D     startvec <- c(Asym = 200, xmid = 725, scal = 350)
##D     nm1 <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal),
##D                   data = Orange,
##D                   fixed = Asym + xmid + scal ~1,
##D                   random = Asym ~1,
##D                   start = startvec)
##D     tidy(nm1)
##D     tidy(nm1, effects = "fixed")
##D     head(augment(nm1, Orange))
##D     glance(nm1)
##D }
## End(Not run)




cleanEx()
nameEx("ordinal_tidiers")
### * ordinal_tidiers

flush(stderr()); flush(stdout())

### Name: tidy.polr
### Title: Tidying methods for ordinal logistic regression models
### Aliases: tidy.polr glance.polr augment.polr ordinal_tidiers tidy.clm
###   tidy.clmm glance.clm glance.clmm augment.clm tidy.svyolr
###   glance.svyolr

### ** Examples

if (require(ordinal)){
  clm_mod <- clm(rating ~ temp * contact, data = wine)
  tidy(clm_mod)
  tidy(clm_mod, conf.int = TRUE)
  tidy(clm_mod, conf.int = TRUE, conf.type = "Wald", exponentiate = TRUE)
  glance(clm_mod)
  augment(clm_mod)

  clm_mod2 <- clm(rating ~ temp, nominal = ~ contact, data = wine)
  tidy(clm_mod2)

  clmm_mod <- clmm(rating ~ temp + contact + (1 | judge), data = wine)
  tidy(clmm_mod)
  glance(clmm_mod)
}
if (require(MASS)) {
  polr_mod <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
  tidy(polr_mod, exponentiate = TRUE, conf.int = TRUE)
  glance(polr_mod)
  augment(polr_mod, type.predict = "class")
}



cleanEx()
nameEx("rowwise_df_tidiers")
### * rowwise_df_tidiers

flush(stderr()); flush(stdout())

### Name: rowwise_df_tidiers
### Title: Tidying methods for rowwise_dfs from dplyr, for tidying each row
###   and recombining the results
### Aliases: rowwise_df_tidiers tidy.rowwise_df tidy_.rowwise_df
###   augment.rowwise_df augment_.rowwise_df glance.rowwise_df
###   glance_.rowwise_df tidy.tbl_df augment.tbl_df glance.tbl_df

### ** Examples


library(dplyr)
regressions <- mtcars %>%
    group_by(cyl) %>%
    do(mod = lm(mpg ~ wt, .))

regressions

regressions %>% tidy(mod)
regressions %>% augment(mod)
regressions %>% glance(mod)

# we can provide additional arguments to the tidying function
regressions %>% tidy(mod, conf.int = TRUE)

# we can also include the original dataset as a "data" argument
# to augment:
regressions <- mtcars %>%
    group_by(cyl) %>%
    do(mod = lm(mpg ~ wt, .), original = (.))

# this allows all the original columns to be included:
regressions %>% augment(mod)  # doesn't include all original
regressions %>% augment(mod, data = original)  # includes all original




cleanEx()
nameEx("rstanarm_tidiers")
### * rstanarm_tidiers

flush(stderr()); flush(stdout())

### Name: rstanarm_tidiers
### Title: Tidying methods for an rstanarm model
### Aliases: rstanarm_tidiers tidy.stanreg glance.stanreg

### ** Examples


## Not run: 
##D fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
##D                   iter = 300, chains = 2)
##D # non-varying ("population") parameters
##D tidy(fit, intervals = TRUE, prob = 0.5)
##D 
##D # hierarchical sd & correlation parameters
##D tidy(fit, parameters = "hierarchical")
##D 
##D # group-specific deviations from "population" parameters
##D tidy(fit, parameters = "varying")
##D 
##D # glance method
##D glance(fit)
##D glance(fit, looic = TRUE, cores = 1)
## End(Not run)




cleanEx()
nameEx("summary_tidiers")
### * summary_tidiers

flush(stderr()); flush(stdout())

### Name: summary_tidiers
### Title: Tidy/glance a(n) summaryDefault object
### Aliases: summary_tidiers tidy.summaryDefault glance.summaryDefault

### ** Examples


v <- rnorm(1000)
s <- summary(v)
s

tidy(s)
glance(s)

v2 <- c(v,NA)
tidy(summary(v2))




cleanEx()
nameEx("tidy.Arima")
### * tidy.Arima

flush(stderr()); flush(stdout())

### Name: tidy.Arima
### Title: Tidy a(n) Arima object
### Aliases: tidy.Arima Arima_tidiers

### ** Examples


fit <- arima(lh, order = c(1, 0, 0))
tidy(fit)
glance(fit)




cleanEx()
nameEx("tidy.Gam")
### * tidy.Gam

flush(stderr()); flush(stdout())

### Name: tidy.Gam
### Title: Tidy a(n) Gam object
### Aliases: tidy.Gam Gam_tidiers

### ** Examples


library(gam)
g <- gam(mpg ~ s(hp, 4) + am + qsec, data = mtcars)
  
tidy(g)
glance(g)




cleanEx()
nameEx("tidy.Kendall")
### * tidy.Kendall

flush(stderr()); flush(stdout())

### Name: tidy.Kendall
### Title: Tidy a(n) Kendall object
### Aliases: tidy.Kendall Kendall_tidiers kendall_tidiers

### ** Examples

library(Kendall)

A <- c(2.5,2.5,2.5,2.5,5,6.5,6.5,10,10,10,10,10,14,14,14,16,17)
B <- c(1,1,1,1,2,1,1,2,1,1,1,1,1,1,2,2,2)

f_res <- Kendall(A, B)
tidy(f_res)

s_res <- MannKendall(B)
tidy(s_res)

t_res <- SeasonalMannKendall(ts(A))
tidy(t_res)




cleanEx()
nameEx("tidy.Mclust")
### * tidy.Mclust

flush(stderr()); flush(stdout())

### Name: tidy.Mclust
### Title: Tidy a(n) Mclust object
### Aliases: tidy.Mclust mclust_tidiers

### ** Examples


library(dplyr) 
library(mclust)
set.seed(27)

centers <- tibble::tibble(
  cluster = factor(1:3), 
  num_points = c(100, 150, 50),  # number points in each cluster
  x1 = c(5, 0, -3),              # x1 coordinate of cluster center
  x2 = c(-1, 1, -2)              # x2 coordinate of cluster center
)

points <- centers %>%
  mutate(
    x1 = purrr::map2(num_points, x1, rnorm),
    x2 = purrr::map2(num_points, x2, rnorm)
  ) %>% 
  select(-num_points, -cluster) %>%
  tidyr::unnest(x1, x2)

m <- mclust::Mclust(points)

tidy(m)
augment(m, points)
glance(m)




cleanEx()
nameEx("tidy.TukeyHSD")
### * tidy.TukeyHSD

flush(stderr()); flush(stdout())

### Name: tidy.TukeyHSD
### Title: Tidy a(n) TukeyHSD object
### Aliases: tidy.TukeyHSD

### ** Examples


fm1 <- aov(breaks ~ wool + tension, data = warpbreaks)
thsd <- TukeyHSD(fm1, "tension", ordered = TRUE)
tidy(thsd)

# may include comparisons on multiple terms
fm2 <- aov(mpg ~ as.factor(gear) * as.factor(cyl), data = mtcars)
tidy(TukeyHSD(fm2))




cleanEx()
nameEx("tidy.aareg")
### * tidy.aareg

flush(stderr()); flush(stdout())

### Name: tidy.aareg
### Title: Tidy a(n) aareg object
### Aliases: tidy.aareg aareg_tidiers

### ** Examples


library(survival)

afit <- aareg(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  dfbeta = TRUE
)

tidy(afit) 




cleanEx()
nameEx("tidy.acf")
### * tidy.acf

flush(stderr()); flush(stdout())

### Name: tidy.acf
### Title: Tidy a(n) acf object
### Aliases: tidy.acf

### ** Examples


tidy(acf(lh, plot = FALSE))
tidy(ccf(mdeaths, fdeaths, plot = FALSE))
tidy(pacf(lh, plot = FALSE))




cleanEx()
nameEx("tidy.anova")
### * tidy.anova

flush(stderr()); flush(stdout())

### Name: tidy.anova
### Title: Tidy a(n) anova object
### Aliases: tidy.anova

### ** Examples


a <- a <- aov(mpg ~ wt + qsec + disp, mtcars)
tidy(a)




cleanEx()
nameEx("tidy.aov")
### * tidy.aov

flush(stderr()); flush(stdout())

### Name: tidy.aov
### Title: Tidy a(n) aov object
### Aliases: tidy.aov

### ** Examples


a <- aov(mpg ~ wt + qsec + disp, mtcars)
tidy(a)




cleanEx()
nameEx("tidy.aovlist")
### * tidy.aovlist

flush(stderr()); flush(stdout())

### Name: tidy.aovlist
### Title: Tidy a(n) aovlist object
### Aliases: tidy.aovlist

### ** Examples


a <- aov(mpg ~ wt + qsec + Error(disp / am), mtcars)
tidy(a)




cleanEx()
nameEx("tidy.betareg")
### * tidy.betareg

flush(stderr()); flush(stdout())

### Name: tidy.betareg
### Title: Tidy a(n) betareg object
### Aliases: tidy.betareg betareg_tidiers

### ** Examples


library(betareg)

data("GasolineYield", package = "betareg")

mod <- betareg(yield ~ batch + temp, data = GasolineYield)

mod
tidy(mod)
tidy(mod, conf.int = TRUE)
tidy(mod, conf.int = TRUE, conf.level = .99)

augment(mod)

glance(mod)




cleanEx()
nameEx("tidy.biglm")
### * tidy.biglm

flush(stderr()); flush(stdout())

### Name: tidy.biglm
### Title: Tidy a(n) biglm object
### Aliases: tidy.biglm

### ** Examples


if (require("biglm", quietly = TRUE)) {
    bfit <- biglm(mpg ~ wt + disp, mtcars)
    tidy(bfit)
    tidy(bfit, conf.int = TRUE)
    tidy(bfit, conf.int = TRUE, conf.level = .9)

    glance(bfit)

    # bigglm: logistic regression
    bgfit <- bigglm(am ~ mpg, mtcars, family = binomial())
    tidy(bgfit)
    tidy(bgfit, exponentiate = TRUE)
    tidy(bgfit, conf.int = TRUE)
    tidy(bgfit, conf.int = TRUE, conf.level = .9)
    tidy(bgfit, conf.int = TRUE, conf.level = .9, exponentiate = TRUE)

    glance(bgfit)
}




cleanEx()
nameEx("tidy.binDesign")
### * tidy.binDesign

flush(stderr()); flush(stdout())

### Name: tidy.binDesign
### Title: Tidy a(n) binDesign object
### Aliases: tidy.binDesign bindesign_tidiers

### ** Examples


if (require("binGroup", quietly = TRUE)) {
    des <- binDesign(nmax = 300, delta = 0.06,
                     p.hyp = 0.1, power = .8)

    glance(des)
    tidy(des)

    # the ggplot2 equivalent of plot(des)
    library(ggplot2)
    ggplot(tidy(des), aes(n, power)) +
        geom_line()
}




cleanEx()
nameEx("tidy.binWidth")
### * tidy.binWidth

flush(stderr()); flush(stdout())

### Name: tidy.binWidth
### Title: Tidy a(n) binWidth object
### Aliases: tidy.binWidth binwidth_tidiers

### ** Examples


if (require("binGroup", quietly = TRUE)) {
    bw <- binWidth(100, .1)
    bw
    tidy(bw)

    library(dplyr)
    d <- expand.grid(n = seq(100, 800, 100),
                     p = .5,
                     method = c("CP", "Blaker", "Score", "Wald"),
                     stringsAsFactors = FALSE) %>%
        group_by(n, p, method) %>%
        do(tidy(binWidth(.$n, .$p, method = .$method)))

    library(ggplot2)
    ggplot(d, aes(n, ci.width, color = method)) +
        geom_line() +
        xlab("Total Observations") +
        ylab("Expected CI Width")
}




cleanEx()
nameEx("tidy.boot")
### * tidy.boot

flush(stderr()); flush(stdout())

### Name: tidy.boot
### Title: Tidy a(n) boot object
### Aliases: tidy.boot boot_tidiers

### ** Examples

if (require("boot")) {
   clotting <- data.frame(
          u = c(5,10,15,20,30,40,60,80,100),
          lot1 = c(118,58,42,35,27,25,21,19,18),
          lot2 = c(69,35,26,21,18,16,13,12,12))

   g1 <- glm(lot2 ~ log(u), data = clotting, family = Gamma)

   bootfun <- function(d, i) {
      coef(update(g1, data= d[i,]))
   }
   bootres <- boot(clotting, bootfun, R = 999)
   tidy(g1, conf.int=TRUE)
   tidy(bootres, conf.int=TRUE)
}




cleanEx()
nameEx("tidy.btergm")
### * tidy.btergm

flush(stderr()); flush(stdout())

### Name: tidy.btergm
### Title: Tidy a(n) btergm object
### Aliases: tidy.btergm btergm_tidiers

### ** Examples


if (require("xergm")) {
    set.seed(1)
    # Using the same simulated example as the xergm package
    # Create 10 random networks with 10 actors
    networks <- list()
    for(i in 1:10){
        mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
        diag(mat) <- 0
        nw <- network::network(mat)
        networks[[i]] <- nw
    }
    # Create 10 matrices as covariates
    covariates <- list()
    for (i in 1:10) {
        mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
        covariates[[i]] <- mat
    }
    # Fit a model where the propensity to form ties depends
    # on the edge covariates, controlling for the number of
    # in-stars
    suppressWarnings(btfit <- btergm(networks ~ edges + istar(2) +
                       edgecov(covariates), R = 100))

    # Show terms, coefficient estimates and errors
    tidy(btfit)

    # Show coefficients as odds ratios with a 99% CI
    tidy(btfit, exponentiate = TRUE, conf.level = 0.99)
}




cleanEx()
nameEx("tidy.cch")
### * tidy.cch

flush(stderr()); flush(stdout())

### Name: tidy.cch
### Title: Tidy a(n) cch object
### Aliases: tidy.cch cch_tidiers

### ** Examples


library(survival)

# examples come from cch documentation
subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel==1|subcoh==1)
ccoh.data <- nwtco[selccoh,]
ccoh.data$subcohort <- subcoh[selccoh]
## central-lab histology
ccoh.data$histol <- factor(ccoh.data$histol,labels=c("FH","UH"))
## tumour stage
ccoh.data$stage <- factor(ccoh.data$stage,labels=c("I","II","III" ,"IV"))
ccoh.data$age <- ccoh.data$age/12 # Age in years

fit.ccP <- cch(Surv(edrel, rel) ~ stage + histol + age, data = ccoh.data,
               subcoh = ~subcohort, id= ~seqno, cohort.size = 4028)

tidy(fit.ccP)

# coefficient plot
library(ggplot2)
ggplot(tidy(fit.ccP), aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0)




cleanEx()
nameEx("tidy.coeftest")
### * tidy.coeftest

flush(stderr()); flush(stdout())

### Name: tidy.coeftest
### Title: Tidy a(n) coeftest object
### Aliases: tidy.coeftest lmtest_tidiers coeftest_tidiers

### ** Examples


if (require("lmtest", quietly = TRUE)) {
    data(Mandible)
    fm <- lm(length ~ age, data=Mandible, subset=(age <= 28))

    lmtest::coeftest(fm)
    tidy(coeftest(fm))
}




cleanEx()
nameEx("tidy.confusionMatrix")
### * tidy.confusionMatrix

flush(stderr()); flush(stdout())

### Name: tidy.confusionMatrix
### Title: Tidy a(n) confusionMatrix object
### Aliases: tidy.confusionMatrix caret_tidiers confusionMatrix_tidiers

### ** Examples


if (requireNamespace("caret", quietly = TRUE)) {

  set.seed(27)
  
  two_class_sample1 <- as.factor(sample(letters[1:2], 100, TRUE))
  two_class_sample2 <- as.factor(sample(letters[1:2], 100, TRUE))
  
  two_class_cm <- caret::confusionMatrix(
    two_class_sample1,
    two_class_sample2
  )
  
  tidy(two_class_cm)
  tidy(two_class_cm, by_class = FALSE)
  
  # multiclass example
  
  six_class_sample1 <- as.factor(sample(letters[1:6], 100, TRUE))
  six_class_sample2 <- as.factor(sample(letters[1:6], 100, TRUE))
  
  six_class_cm <- caret::confusionMatrix(
    six_class_sample1,
    six_class_sample2
  )
  
  tidy(six_class_cm)
  tidy(six_class_cm, by_class = FALSE)
}




cleanEx()
nameEx("tidy.coxph")
### * tidy.coxph

flush(stderr()); flush(stdout())

### Name: tidy.coxph
### Title: Tidy a(n) coxph object
### Aliases: tidy.coxph coxph_tidiers

### ** Examples


library(survival)

cfit <- coxph(Surv(time, status) ~ age + sex, lung)

tidy(cfit)
tidy(cfit, exponentiate = TRUE)

lp <- augment(cfit, lung)
risks <- augment(cfit, lung, type.predict = "risk")
expected <- augment(cfit, lung, type.predict = "expected")

glance(cfit)

# also works on clogit models
resp <- levels(logan$occupation)
n <- nrow(logan)
indx <- rep(1:n, length(resp))
logan2 <- data.frame(
  logan[indx,],
  id = indx,
  tocc = factor(rep(resp, each=n))
)

logan2$case <- (logan2$occupation == logan2$tocc)

cl <- clogit(case ~ tocc + tocc:education + strata(id), logan2)
tidy(cl)
glance(cl)

library(ggplot2)

ggplot(lp, aes(age, .fitted, color = sex)) +
  geom_point()

ggplot(risks, aes(age, .fitted, color = sex)) + 
  geom_point()
  
ggplot(expected, aes(time, .fitted, color = sex)) + 
  geom_point()





cleanEx()
nameEx("tidy.cv.glmnet")
### * tidy.cv.glmnet

flush(stderr()); flush(stdout())

### Name: tidy.cv.glmnet
### Title: Tidy a(n) cv.glmnet object
### Aliases: tidy.cv.glmnet

### ** Examples


if (requireNamespace("glmnet", quietly = TRUE)) {

    library(glmnet)
    set.seed(27)

    nobs <- 100
    nvar <- 50
    real <- 5

    x <- matrix(rnorm(nobs * nvar), nobs, nvar)
    beta <- c(rnorm(real, 0, 1), rep(0, nvar - real))
    y <- c(t(beta) %*% t(x)) + rnorm(nvar, sd = 3)

    cvfit1 <- cv.glmnet(x,y)

    tidy(cvfit1)
    glance(cvfit1)

    library(ggplot2)
    tidied_cv <- tidy(cvfit1)
    glance_cv <- glance(cvfit1)

    # plot of MSE as a function of lambda
    g <- ggplot(tidied_cv, aes(lambda, estimate)) + geom_line() + scale_x_log10()
    g

    # plot of MSE as a function of lambda with confidence ribbon
    g <- g + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .25)
    g

    # plot of MSE as a function of lambda with confidence ribbon and choices
    # of minimum lambda marked
    g <- g + geom_vline(xintercept = glance_cv$lambda.min) +
        geom_vline(xintercept = glance_cv$lambda.1se, lty = 2)
    g

    # plot of number of zeros for each choice of lambda
    ggplot(tidied_cv, aes(lambda, nzero)) + geom_line() + scale_x_log10()

    # coefficient plot with min lambda shown
    tidied <- tidy(cvfit1$glmnet.fit)
    ggplot(tidied, aes(lambda, estimate, group = term)) + scale_x_log10() +
        geom_line() +
        geom_vline(xintercept = glance_cv$lambda.min) +
        geom_vline(xintercept = glance_cv$lambda.1se, lty = 2)
}




cleanEx()
nameEx("tidy.dist")
### * tidy.dist

flush(stderr()); flush(stdout())

### Name: tidy.dist
### Title: Tidy a(n) dist object
### Aliases: tidy.dist

### ** Examples


iris_dist <- dist(t(iris[, 1:4]))
iris_dist

tidy(iris_dist)
tidy(iris_dist, upper = TRUE)
tidy(iris_dist, diagonal = TRUE)




cleanEx()
nameEx("tidy.ergm")
### * tidy.ergm

flush(stderr()); flush(stdout())

### Name: tidy.ergm
### Title: Tidy a(n) ergm object
### Aliases: tidy.ergm ergm_tidiers

### ** Examples


library(ergm)
# Using the same example as the ergm package
# Load the Florentine marriage network data
data(florentine)

# Fit a model where the propensity to form ties between
# families depends on the absolute difference in wealth
gest <- ergm(flomarriage ~ edges + absdiff("wealth"))

# Show terms, coefficient estimates and errors
tidy(gest)

# Show coefficients as odds ratios with a 99% CI
tidy(gest, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.99)

# Take a look at likelihood measures and other
# control parameters used during MCMC estimation
glance(gest)
glance(gest, deviance = TRUE)
glance(gest, mcmc = TRUE)




cleanEx()
nameEx("tidy.factanal")
### * tidy.factanal

flush(stderr()); flush(stdout())

### Name: tidy.factanal
### Title: Tidy a(n) factanal object
### Aliases: tidy.factanal factanal_tidiers

### ** Examples


mod <- factanal(mtcars, 3, scores = "regression")

glance(mod)
tidy(mod)
augment(mod)
augment(mod, mtcars)




cleanEx()
nameEx("tidy.felm")
### * tidy.felm

flush(stderr()); flush(stdout())

### Name: tidy.felm
### Title: Tidy a(n) felm object
### Aliases: tidy.felm felm_tidiers lfe_tidiers

### ** Examples


if (require("lfe", quietly = TRUE)) {

    library(lfe)
    
    N=1e2
    DT <- data.frame(
      id = sample(5, N, TRUE),
      v1 =  sample(5, N, TRUE),
      v2 =  sample(1e6, N, TRUE),
      v3 =  sample(round(runif(100,max=100),4), N, TRUE),
      v4 =  sample(round(runif(100,max=100),4), N, TRUE)
    )

    result_felm <- felm(v2~v3, DT)
    tidy(result_felm)
    augment(result_felm)
    result_felm <- felm(v2~v3|id+v1, DT)
    tidy(result_felm, fe = TRUE)
    augment(result_felm)
    v1<-DT$v1
    v2 <- DT$v2
    v3 <- DT$v3
    id <- DT$id
    result_felm <- felm(v2~v3|id+v1)
    tidy(result_felm)
    augment(result_felm)
    glance(result_felm)
}




cleanEx()
nameEx("tidy.fitdistr")
### * tidy.fitdistr

flush(stderr()); flush(stdout())

### Name: tidy.fitdistr
### Title: Tidy a(n) fitdistr object
### Aliases: tidy.fitdistr fitdistr_tidiers

### ** Examples


set.seed(2015)
x <- rnorm(100, 5, 2)

library(MASS)
fit <- fitdistr(x, dnorm, list(mean = 3, sd = 1))

tidy(fit)
glance(fit)




cleanEx()
nameEx("tidy.ftable")
### * tidy.ftable

flush(stderr()); flush(stdout())

### Name: tidy.ftable
### Title: Tidy a(n) ftable object
### Aliases: tidy.ftable

### ** Examples


tidy(ftable(Titanic, row.vars = 1:3))




cleanEx()
nameEx("tidy.gamlss")
### * tidy.gamlss

flush(stderr()); flush(stdout())

### Name: tidy.gamlss
### Title: Tidy a(n) gamlss object
### Aliases: tidy.gamlss

### ** Examples


library(gamlss)

g <- gamlss(
  y ~ pb(x),
  sigma.fo = ~ pb(x),
  family = BCT,
  data = abdom,
  method = mixed(1, 20)
)

tidy(g)




cleanEx()
nameEx("tidy.garch")
### * tidy.garch

flush(stderr()); flush(stdout())

### Name: tidy.garch
### Title: Tidy a(n) garch object
### Aliases: tidy.garch garch_tidiers

### ** Examples


library(tseries)

data(EuStockMarkets)
dax <- diff(log(EuStockMarkets))[,"DAX"]
dax.garch <- garch(dax)
dax.garch

tidy(dax.garch)
glance(dax.garch)




cleanEx()
nameEx("tidy.glht")
### * tidy.glht

flush(stderr()); flush(stdout())

### Name: tidy.glht
### Title: Tidy a(n) glht object
### Aliases: tidy.glht multcomp_tidiers

### ** Examples


if (require("multcomp") && require("ggplot2")) {

    library(multcomp)
    library(ggplot2)
    
    amod <- aov(breaks ~ wool + tension, data = warpbreaks)
    wht <- glht(amod, linfct = mcp(tension = "Tukey"))

    tidy(wht)
    ggplot(wht, aes(lhs, estimate)) + geom_point()

    CI <- confint(wht)
    tidy(CI)
    ggplot(CI, aes(lhs, estimate, ymin = lwr, ymax = upr)) +
       geom_pointrange()

    tidy(summary(wht))
    ggplot(mapping = aes(lhs, estimate)) +
       geom_linerange(aes(ymin = lwr, ymax = upr), data = CI) +
       geom_point(aes(size = p), data = summary(wht)) +
       scale_size(trans = "reverse")

    cld <- cld(wht)
    tidy(cld)
}




cleanEx()
nameEx("tidy.glmRob")
### * tidy.glmRob

flush(stderr()); flush(stdout())

### Name: tidy.glmRob
### Title: Tidy a(n) glmRob object
### Aliases: tidy.glmRob

### ** Examples


library(robust)
m <- lmRob(mpg ~ wt, data = mtcars)

tidy(m)
augment(m)
glance(m)

gm <- glmRob(am ~ wt, data = mtcars, family = "binomial")
glance(gm)




cleanEx()
nameEx("tidy.glmnet")
### * tidy.glmnet

flush(stderr()); flush(stdout())

### Name: tidy.glmnet
### Title: Tidy a(n) glmnet object
### Aliases: tidy.glmnet glmnet_tidiers

### ** Examples


if (requireNamespace("glmnet", quietly = TRUE)) {

    library(glmnet)
    
    set.seed(2014)
    x <- matrix(rnorm(100*20),100,20)
    y <- rnorm(100)
    fit1 <- glmnet(x,y)

    tidy(fit1)
    glance(fit1)

    library(dplyr)
    library(ggplot2)

    tidied <- tidy(fit1) %>% filter(term != "(Intercept)")

    ggplot(tidied, aes(step, estimate, group = term)) + geom_line()
    ggplot(tidied, aes(lambda, estimate, group = term)) +
        geom_line() + scale_x_log10()

    ggplot(tidied, aes(lambda, dev.ratio)) + geom_line()

    # works for other types of regressions as well, such as logistic
    g2 <- sample(1:2, 100, replace=TRUE)
    fit2 <- glmnet(x, g2, family="binomial")
    tidy(fit2)
}




cleanEx()
nameEx("tidy.gmm")
### * tidy.gmm

flush(stderr()); flush(stdout())

### Name: tidy.gmm
### Title: Tidy a(n) gmm object
### Aliases: tidy.gmm gmm_tidiers

### ** Examples


if (requireNamespace("gmm", quietly = TRUE)) {

  library(gmm)
  
  # examples come from the "gmm" package
  ## CAPM test with GMM
  data(Finance)
  r <- Finance[1:300, 1:10]
  rm <- Finance[1:300, "rm"]
  rf <- Finance[1:300, "rf"]

  z <- as.matrix(r-rf)
  t <- nrow(z)
  zm <- rm-rf
  h <- matrix(zm, t, 1)
  res <- gmm(z ~ zm, x = h)

  # tidy result
  tidy(res)
  tidy(res, conf.int = TRUE)
  tidy(res, conf.int = TRUE, conf.level = .99)

  # coefficient plot
  library(ggplot2)
  library(dplyr)
  tidy(res, conf.int = TRUE) %>%
    mutate(variable = reorder(variable, estimate)) %>%
    ggplot(aes(estimate, variable)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
    facet_wrap(~ term) +
    geom_vline(xintercept = 0, color = "red", lty = 2)

  # from a function instead of a matrix
  g <- function(theta, x) {
  	e <- x[,2:11] - theta[1] - (x[,1] - theta[1]) %*% matrix(theta[2:11], 1, 10)
  	gmat <- cbind(e, e*c(x[,1]))
  	return(gmat) }

  x <- as.matrix(cbind(rm, r))
  res_black <- gmm(g, x = x, t0 = rep(0, 11))

  tidy(res_black)
  tidy(res_black, conf.int = TRUE)

  ## APT test with Fama-French factors and GMM

  f1 <- zm
  f2 <- Finance[1:300, "hml"] - rf
  f3 <- Finance[1:300, "smb"] - rf
  h <- cbind(f1, f2, f3)
  res2 <- gmm(z ~ f1 + f2 + f3, x = h)

  td2 <- tidy(res2, conf.int = TRUE)
  td2

  # coefficient plot
  td2 %>%
    mutate(variable = reorder(variable, estimate)) %>%
    ggplot(aes(estimate, variable)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
    facet_wrap(~ term) +
    geom_vline(xintercept = 0, color = "red", lty = 2)
}




cleanEx()
nameEx("tidy.htest")
### * tidy.htest

flush(stderr()); flush(stdout())

### Name: tidy.htest
### Title: Tidy/glance a(n) htest object
### Aliases: tidy.htest htest_tidiers glance.htest

### ** Examples


tt <- t.test(rnorm(10))
tidy(tt)
glance(tt)  # same output for all htests

tt <- t.test(mpg ~ am, data = mtcars)
tidy(tt)

wt <- wilcox.test(mpg ~ am, data = mtcars, conf.int = TRUE, exact = FALSE)
tidy(wt)

ct <- cor.test(mtcars$wt, mtcars$mpg)
tidy(ct)

chit <- chisq.test(xtabs(Freq ~ Sex + Class, data = as.data.frame(Titanic)))
tidy(chit)
augment(chit)




cleanEx()
nameEx("tidy.ivreg")
### * tidy.ivreg

flush(stderr()); flush(stdout())

### Name: tidy.ivreg
### Title: Tidy a(n) ivreg object
### Aliases: tidy.ivreg ivreg_tidiers aer_tidiers

### ** Examples


library(AER)

data("CigarettesSW", package = "AER")
ivr <- ivreg(
  log(packs) ~ income | population,
  data = CigarettesSW,
  subset = year == "1995"
)

summary(ivr)

tidy(ivr)
tidy(ivr, conf.int = TRUE)
tidy(ivr, conf.int = TRUE, exponentiate = TRUE)

augment(ivr)

glance(ivr)




cleanEx()
nameEx("tidy.kappa")
### * tidy.kappa

flush(stderr()); flush(stdout())

### Name: tidy.kappa
### Title: Tidy a(n) kappa object
### Aliases: tidy.kappa kappa_tidiers psych_tidiers

### ** Examples


library(psych)

rater1 = 1:9
rater2 = c(1, 3, 1, 6, 1, 5, 5, 6, 7)
ck <- cohen.kappa(cbind(rater1, rater2))

tidy(ck)

# graph the confidence intervals
library(ggplot2)
ggplot(tidy(ck), aes(estimate, type)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high))




cleanEx()
nameEx("tidy.kde")
### * tidy.kde

flush(stderr()); flush(stdout())

### Name: tidy.kde
### Title: Tidy a(n) kde object
### Aliases: tidy.kde kde_tidiers ks_tidiers

### ** Examples


if (requireNamespace("ks", quietly = TRUE)) {
  
  library(ks)
  
  dat <- replicate(2, rnorm(100))
  k <- kde(dat)

  td <- tidy(k)
  td

  library(ggplot2)
  ggplot(td, aes(x1, x2, fill = estimate)) +
    geom_tile() +
    theme_void()

  # also works with 3 dimensions
  dat3 <- replicate(3, rnorm(100))
  k3 <- kde(dat3)

  td3 <- tidy(k3)
  td3
}




cleanEx()
nameEx("tidy.lavaan")
### * tidy.lavaan

flush(stderr()); flush(stdout())

### Name: tidy.lavaan
### Title: Tidy a(n) lavaan object
### Aliases: tidy.lavaan lavaan_tidiers sem_tidiers cfa_tidiers

### ** Examples


if (require("lavaan")) {

 library(lavaan)
 
 cfa.fit <- cfa('F =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9',
                data = HolzingerSwineford1939, group = "school")
 tidy(cfa.fit)
}




cleanEx()
nameEx("tidy.lm")
### * tidy.lm

flush(stderr()); flush(stdout())

### Name: tidy.lm
### Title: Tidy a(n) lm object
### Aliases: tidy.lm lm_tidiers tidy.summary.lm

### ** Examples


library(ggplot2)
library(dplyr)

mod <- lm(mpg ~ wt + qsec, data = mtcars)

tidy(mod)
glance(mod)

# coefficient plot
d <- tidy(mod) %>% 
  mutate(
    low = estimate - std.error,
    high = estimate + std.error
  )
  
ggplot(d, aes(estimate, term, xmin = low, xmax = high, height = 0)) +
     geom_point() +
     geom_vline(xintercept = 0) +
     geom_errorbarh()

augment(mod)
augment(mod, mtcars)

# predict on new data
newdata <- mtcars %>% head(6) %>% mutate(wt = wt + 1)
augment(mod, newdata = newdata)

au <- augment(mod, data = mtcars)

ggplot(au, aes(.hat, .std.resid)) +
  geom_vline(size = 2, colour = "white", xintercept = 0) +
  geom_hline(size = 2, colour = "white", yintercept = 0) +
  geom_point() + geom_smooth(se = FALSE)

plot(mod, which = 6)
ggplot(au, aes(.hat, .cooksd)) +
  geom_vline(xintercept = 0, colour = NA) +
  geom_abline(slope = seq(0, 3, by = 0.5), colour = "white") +
  geom_smooth(se = FALSE) +
  geom_point()

# column-wise models
a <- matrix(rnorm(20), nrow = 10)
b <- a + rnorm(length(a))
result <- lm(b ~ a)
tidy(result)




cleanEx()
nameEx("tidy.lmRob")
### * tidy.lmRob

flush(stderr()); flush(stdout())

### Name: tidy.lmRob
### Title: Tidy a(n) lmRob object
### Aliases: tidy.lmRob robust_tidiers

### ** Examples


library(robust)
m <- lmRob(mpg ~ wt, data = mtcars)

tidy(m)
augment(m)
glance(m)

gm <- glmRob(am ~ wt, data = mtcars, family = "binomial")
glance(gm)




cleanEx()
nameEx("tidy.lmodel2")
### * tidy.lmodel2

flush(stderr()); flush(stdout())

### Name: tidy.lmodel2
### Title: Tidy a(n) lmodel2 object
### Aliases: tidy.lmodel2 lmodel2_tidiers

### ** Examples


if (require("lmodel2", quietly = TRUE)) {

  library(lmodel2)
  
  data(mod2ex2)
  Ex2.res <- lmodel2(Prey ~ Predators, data=mod2ex2, "relative", "relative", 99)
  Ex2.res

  tidy(Ex2.res)
  glance(Ex2.res)

  # this allows coefficient plots with ggplot2
  library(ggplot2)
  ggplot(tidy(Ex2.res), aes(estimate, term, color = method)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high))
}




cleanEx()
nameEx("tidy.manova")
### * tidy.manova

flush(stderr()); flush(stdout())

### Name: tidy.manova
### Title: Tidy a(n) manova object
### Aliases: tidy.manova

### ** Examples


npk2 <- within(npk, foo <- rnorm(24))
m <- manova(cbind(yield, foo) ~ block + N * P * K, npk2)
tidy(m) 




cleanEx()
nameEx("tidy.map")
### * tidy.map

flush(stderr()); flush(stdout())

### Name: tidy.map
### Title: Tidy a(n) map object
### Aliases: tidy.map maps_tidiers

### ** Examples


if (require("maps") && require("ggplot2")) {
    
    library(maps)
    library(ggplot2)

    ca <- map("county", "ca", plot = FALSE, fill = TRUE)
    tidy(ca)
    qplot(long, lat, data = ca, geom = "polygon", group = group)

    tx <- map("county", "texas", plot = FALSE, fill = TRUE)
    tidy(tx)
    qplot(long, lat, data = tx, geom = "polygon", group = group,
          colour = I("white"))
}




cleanEx()
nameEx("tidy.mjoint")
### * tidy.mjoint

flush(stderr()); flush(stdout())

### Name: tidy.mjoint
### Title: Tidy a(n) mjoint object
### Aliases: tidy.mjoint mjoint_tidiers joinerml_tidiers

### ** Examples

## Not run: 
##D # Fit a joint model with bivariate longitudinal outcomes
##D library(joineRML)
##D data(heart.valve)
##D hvd <- heart.valve[!is.na(heart.valve$log.grad) &
##D                        !is.na(heart.valve$log.lvmi) &
##D                        heart.valve$num <= 50, ]
##D fit <- mjoint(
##D     formLongFixed = list(
##D         "grad" = log.grad ~ time + sex + hs,
##D         "lvmi" = log.lvmi ~ time + sex
##D     ),
##D     formLongRandom = list(
##D         "grad" = ~ 1 | num,
##D         "lvmi" = ~ time | num
##D     ),
##D     formSurv = Surv(fuyrs, status) ~ age,
##D     data = hvd,
##D     inits = list("gamma" = c(0.11, 1.51, 0.80)),
##D     timeVar = "time"
##D )
##D 
##D # Extract the survival fixed effects
##D tidy(fit)
##D 
##D # Extract the longitudinal fixed effects
##D tidy(fit, component = "longitudinal")
##D 
##D # Extract the survival fixed effects with confidence intervals
##D tidy(fit, ci = TRUE)
##D 
##D # Extract the survival fixed effects with confidence intervals based
##D # on bootstrapped standard errors
##D bSE <- bootSE(fit, nboot = 5, safe.boot = TRUE)
##D tidy(fit, boot_se = bSE, ci = TRUE)
##D 
##D # Augment original data with fitted longitudinal values and residuals
##D hvd2 <- augment(fit)
##D 
##D # Extract model statistics
##D glance(fit)
## End(Not run)




cleanEx()
nameEx("tidy.mle2")
### * tidy.mle2

flush(stderr()); flush(stdout())

### Name: tidy.mle2
### Title: Tidy a(n) mle2 object
### Aliases: tidy.mle2 mle2_tidiers bbmle_tidiers

### ** Examples


if (require("bbmle", quietly = TRUE)) {
  x <- 0:10
  y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
  d <- data.frame(x,y)

  fit <- mle2(y ~ dpois(lambda = ymean),
              start = list(ymean = mean(y)), data = d)

  tidy(fit)
}




cleanEx()
nameEx("tidy.muhaz")
### * tidy.muhaz

flush(stderr()); flush(stdout())

### Name: tidy.muhaz
### Title: Tidy a(n) muhaz object
### Aliases: tidy.muhaz muhaz_tidiers

### ** Examples

if (require("muhaz", quietly = TRUE)) {
  data(ovarian, package="survival")
  x <- muhaz::muhaz(ovarian$futime, ovarian$fustat)
  tidy(x)
  glance(x)
}




cleanEx()
nameEx("tidy.multinom")
### * tidy.multinom

flush(stderr()); flush(stdout())

### Name: tidy.multinom
### Title: Tidying methods for multinomial logistic regression models
### Aliases: tidy.multinom multinom_tidiers nnet_tidiers

### ** Examples


if (require(nnet) & require(MASS)){
  library(nnet)
  library(MASS)
  
  example(birthwt)
  bwt.mu <- multinom(low ~ ., bwt)
  tidy(bwt.mu)
  glance(bwt.mu)

  #* This model is a truly terrible model
  #* but it should show you what the output looks
  #* like in a multinomial logistic regression

  fit.gear <- multinom(gear ~ mpg + factor(am), data = mtcars)
  tidy(fit.gear)
  glance(fit.gear)
}




cleanEx()
nameEx("tidy.nls")
### * tidy.nls

flush(stderr()); flush(stdout())

### Name: tidy.nls
### Title: Tidy a(n) nls object
### Aliases: tidy.nls nls_tidiers

### ** Examples


n <- nls(mpg ~ k * e ^ wt, data = mtcars, start = list(k = 1, e = 2))

tidy(n)
augment(n)
glance(n)

library(ggplot2)
ggplot(augment(n), aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = .fitted))

newdata <- head(mtcars)
newdata$wt <- newdata$wt + 1
augment(n, newdata = newdata)




cleanEx()
nameEx("tidy.orcutt")
### * tidy.orcutt

flush(stderr()); flush(stdout())

### Name: tidy.orcutt
### Title: Tidy a(n) orcutt object
### Aliases: tidy.orcutt orcutt_tidiers

### ** Examples


reg <- lm(mpg ~ wt + qsec + disp, mtcars)
tidy(reg)

if (require("orcutt", quietly = TRUE)) {
  co <- cochrane.orcutt(reg)
  co

  tidy(co)
  glance(co)
}




cleanEx()
nameEx("tidy.pairwise.htest")
### * tidy.pairwise.htest

flush(stderr()); flush(stdout())

### Name: tidy.pairwise.htest
### Title: Tidy a(n) pairwise.htest object
### Aliases: tidy.pairwise.htest

### ** Examples


attach(airquality)
Month <- factor(Month, labels = month.abb[5:9])
ptt <- pairwise.t.test(Ozone, Month)
tidy(ptt)

attach(iris)
ptt2 <- pairwise.t.test(Petal.Length, Species)
tidy(ptt2)

tidy(pairwise.t.test(Petal.Length, Species, alternative = "greater"))
tidy(pairwise.t.test(Petal.Length, Species, alternative = "less"))

tidy(pairwise.wilcox.test(Petal.Length, Species))




cleanEx()
nameEx("tidy.plm")
### * tidy.plm

flush(stderr()); flush(stdout())

### Name: tidy.plm
### Title: Tidy a(n) plm object
### Aliases: tidy.plm plm_tidiers

### ** Examples


library(plm)

data("Produc", package = "plm")
zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
          data = Produc, index = c("state","year"))

summary(zz)

tidy(zz)
tidy(zz, conf.int = TRUE)
tidy(zz, conf.int = TRUE, conf.level = .9)

augment(zz)
glance(zz)




cleanEx()
nameEx("tidy.poLCA")
### * tidy.poLCA

flush(stderr()); flush(stdout())

### Name: tidy.poLCA
### Title: Tidy a(n) poLCA object
### Aliases: tidy.poLCA poLCA_tidiers

### ** Examples


if (require("poLCA", quietly = TRUE)) {
  library(poLCA)
  library(dplyr)

  data(values)
  f <- cbind(A, B, C, D)~1
  M1 <- poLCA(f, values, nclass = 2, verbose = FALSE)

  M1
  tidy(M1)
  augment(M1)
  glance(M1)

  library(ggplot2)

  ggplot(tidy(M1), aes(factor(class), estimate, fill = factor(outcome))) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(~ variable)

  set.seed(2016)
  # compare multiple
  mods <- tibble(nclass = 1:3) %>%
    group_by(nclass) %>%
    do(mod = poLCA(f, values, nclass = .$nclass, verbose = FALSE))

  # compare log-likelihood and/or AIC, BIC
  mods %>%
    glance(mod)

  ## Three-class model with a single covariate.

  data(election)
  f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
               MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY
  nes2a <- poLCA(f2a, election, nclass = 3, nrep = 5, verbose = FALSE)

  td <- tidy(nes2a)
  td

  # show

  ggplot(td, aes(outcome, estimate, color = factor(class), group = class)) +
    geom_line() +
    facet_wrap(~ variable, nrow = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  au <- augment(nes2a)
  au
  au %>%
    count(.class)

  # if the original data is provided, it leads to NAs in new columns
  # for rows that weren't predicted
  au2 <- augment(nes2a, data = election)
  au2
  dim(au2)
}




cleanEx()
nameEx("tidy.power.htest")
### * tidy.power.htest

flush(stderr()); flush(stdout())

### Name: tidy.power.htest
### Title: Tidy a(n) power.htest object
### Aliases: tidy.power.htest

### ** Examples


ptt <- power.t.test(n = 2:30, delta = 1)
tidy(ptt)

library(ggplot2)

ggplot(tidy(ptt), aes(n, power)) +
  geom_line()




cleanEx()
nameEx("tidy.prcomp")
### * tidy.prcomp

flush(stderr()); flush(stdout())

### Name: tidy.prcomp
### Title: Tidy a(n) prcomp object
### Aliases: tidy.prcomp prcomp_tidiers

### ** Examples


pc <- prcomp(USArrests, scale = TRUE)

# information about rotation
tidy(pc)

# information about samples (states)
tidy(pc, "samples")

# information about PCs
tidy(pc, "pcs")

# state map
library(dplyr)
library(ggplot2)

pc %>%
  tidy(matrix = "samples") %>%
  mutate(region = tolower(row)) %>%
  inner_join(map_data("state"), by = "region") %>%
  ggplot(aes(long, lat, group = group, fill = value)) +
  geom_polygon() +
  facet_wrap(~ PC) +
  theme_void() +
  ggtitle("Principal components of arrest data")

au <- augment(pc, data = USArrests)
au

ggplot(au, aes(.fittedPC1, .fittedPC2)) +
  geom_point() +
  geom_text(aes(label = .rownames), vjust = 1, hjust = 1)




cleanEx()
nameEx("tidy.pyears")
### * tidy.pyears

flush(stderr()); flush(stdout())

### Name: tidy.pyears
### Title: Tidy a(n) pyears object
### Aliases: tidy.pyears pyears_tidiers

### ** Examples


library(survival)

temp.yr  <- tcut(mgus$dxyr, 55:92, labels=as.character(55:91))
temp.age <- tcut(mgus$age, 34:101, labels=as.character(34:100))
ptime <- ifelse(is.na(mgus$pctime), mgus$futime, mgus$pctime)
pstat <- ifelse(is.na(mgus$pctime), 0, 1)
pfit <- pyears(Surv(ptime/365.25, pstat) ~ temp.yr + temp.age + sex,  mgus,
               data.frame=TRUE)
tidy(pfit)
glance(pfit)

# if data.frame argument is not given, different information is present in
# output
pfit2 <- pyears(Surv(ptime/365.25, pstat) ~ temp.yr + temp.age + sex,  mgus)
tidy(pfit2)
glance(pfit2)




cleanEx()
nameEx("tidy.rcorr")
### * tidy.rcorr

flush(stderr()); flush(stdout())

### Name: tidy.rcorr
### Title: Tidy a(n) rcorr object
### Aliases: tidy.rcorr rcorr_tidiers Hmisc_tidiers

### ** Examples


if (requireNamespace("Hmisc", quietly = TRUE)) {

    library(Hmisc)
    
    mat <- replicate(52, rnorm(100))
    # add some NAs
    mat[sample(length(mat), 2000)] <- NA
    # also column names
    colnames(mat) <- c(LETTERS, letters)

    rc <- rcorr(mat)

    td <- tidy(rc)
    td

    library(ggplot2)
    ggplot(td, aes(p.value)) +
        geom_histogram(binwidth = .1)

    ggplot(td, aes(estimate, p.value)) +
        geom_point() +
        scale_y_log10()
}




cleanEx()
nameEx("tidy.ridgelm")
### * tidy.ridgelm

flush(stderr()); flush(stdout())

### Name: tidy.ridgelm
### Title: Tidy a(n) ridgelm object
### Aliases: tidy.ridgelm ridgelm_tidiers

### ** Examples


names(longley)[1] <- "y"
fit1 <- MASS::lm.ridge(y ~ ., longley)
tidy(fit1)

fit2 <- MASS::lm.ridge(y ~ ., longley, lambda = seq(0.001, .05, .001))
td2 <- tidy(fit2)
g2 <- glance(fit2)

# coefficient plot
library(ggplot2)
ggplot(td2, aes(lambda, estimate, color = term)) +
  geom_line()

# GCV plot
ggplot(td2, aes(lambda, GCV)) +
  geom_line()

# add line for the GCV minimizing estimate
ggplot(td2, aes(lambda, GCV)) + 
  geom_line() +
  geom_vline(xintercept = g2$lambdaGCV, col = "red", lty = 2)




cleanEx()
nameEx("tidy.roc")
### * tidy.roc

flush(stderr()); flush(stdout())

### Name: tidy.roc
### Title: Tidy a(n) roc object
### Aliases: tidy.roc auc_tidiers roc_tidiers

### ** Examples


if (require("AUC", quietly = TRUE)) {
  data(churn)
  r <- roc(churn$predictions,churn$labels)

  td <- tidy(r)
  td

  library(ggplot2)
  
  ggplot(td, aes(fpr, tpr)) +
    geom_line()

  # compare the ROC curves for two prediction algorithms
  
  library(dplyr)
  library(tidyr)

  rocs <- churn %>%
    gather(algorithm, value, -labels) %>%
    nest(-algorithm) %>% 
    mutate(tidy_roc = purrr::map(data, ~tidy(roc(.x$value, .x$labels)))) %>% 
    unnest(tidy_roc)

  ggplot(rocs, aes(fpr, tpr, color = algorithm)) +
    geom_line()
}




cleanEx()
nameEx("tidy.spec")
### * tidy.spec

flush(stderr()); flush(stdout())

### Name: tidy.spec
### Title: Tidy a(n) spec object
### Aliases: tidy.spec

### ** Examples


spc <- spectrum(lh)
tidy(spc)

library(ggplot2)
ggplot(tidy(spc), aes(freq, spec)) +
  geom_line()




cleanEx()
nameEx("tidy.speedlm")
### * tidy.speedlm

flush(stderr()); flush(stdout())

### Name: tidy.speedlm
### Title: Tidy a(n) speedlm object
### Aliases: tidy.speedlm speedlm_tidiers speedglm_tidiers

### ** Examples


mod <- speedglm::speedlm(mpg ~ wt + qsec, data = mtcars)

tidy(mod)
glance(mod)
augment(mod)




cleanEx()
nameEx("tidy.survdiff")
### * tidy.survdiff

flush(stderr()); flush(stdout())

### Name: tidy.survdiff
### Title: Tidy a(n) survdiff object
### Aliases: tidy.survdiff survdiff_tidiers

### ** Examples


library(survival)

s <- survdiff(
  Surv(time, status) ~ pat.karno + strata(inst),
  data = lung
)

tidy(s)
glance(s)




cleanEx()
nameEx("tidy.survexp")
### * tidy.survexp

flush(stderr()); flush(stdout())

### Name: tidy.survexp
### Title: Tidy a(n) survexp object
### Aliases: tidy.survexp sexpfit_tidiers survexp_tidiers

### ** Examples


library(survival)
sexpfit <- survexp(
  futime ~ 1,
  rmap = list(
    sex = "male",
    year = accept.dt,
    age = (accept.dt - birth.dt)
  ),
  method = 'conditional',
  data = jasa
)

tidy(sexpfit)
glance(sexpfit)




cleanEx()
nameEx("tidy.survfit")
### * tidy.survfit

flush(stderr()); flush(stdout())

### Name: tidy.survfit
### Title: Tidy a(n) survfit object
### Aliases: tidy.survfit survfit_tidiers

### ** Examples


library(survival)
cfit <- coxph(Surv(time, status) ~ age + sex, lung)
sfit <- survfit(cfit)

tidy(sfit)
glance(sfit)

library(ggplot2)
ggplot(tidy(sfit), aes(time, estimate)) + geom_line() +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25)

# multi-state
fitCI <- survfit(Surv(stop, status * as.numeric(event), type = "mstate") ~ 1,
              data = mgus1, subset = (start == 0))
td_multi <- tidy(fitCI)
td_multi

ggplot(td_multi, aes(time, estimate, group = state)) +
    geom_line(aes(color = state)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .25)




cleanEx()
nameEx("tidy.survreg")
### * tidy.survreg

flush(stderr()); flush(stdout())

### Name: tidy.survreg
### Title: Tidy a(n) survreg object
### Aliases: tidy.survreg survreg_tidiers

### ** Examples


library(survival)

sr <- survreg(
  Surv(futime, fustat) ~ ecog.ps + rx,
  ovarian,
  dist = "exponential"
)

td <- tidy(sr)
augment(sr, ovarian)
glance(sr)

# coefficient plot
library(ggplot2)
ggplot(td, aes(estimate, term)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0)




cleanEx()
nameEx("tidy.table")
### * tidy.table

flush(stderr()); flush(stdout())

### Name: tidy.table
### Title: Tidy a(n) table object
### Aliases: tidy.table

### ** Examples


tab <- with(airquality, table(Temp = cut(Temp, quantile(Temp)), Month))
tidy(tab)




cleanEx()
nameEx("tidy.ts")
### * tidy.ts

flush(stderr()); flush(stdout())

### Name: tidy.ts
### Title: Tidy a(n) ts object
### Aliases: tidy.ts

### ** Examples


set.seed(678)

tidy(ts(1:10, frequency = 4, start = c(1959, 2)))

z <- ts(matrix(rnorm(300), 100, 3), start = c(1961, 1), frequency = 12)
colnames(z) <- c("Aa", "Bb", "Cc")
tidy(z)




cleanEx()
nameEx("tidy.zoo")
### * tidy.zoo

flush(stderr()); flush(stdout())

### Name: tidy.zoo
### Title: Tidy a(n) zoo object
### Aliases: tidy.zoo zoo_tidiers

### ** Examples


library(zoo)
library(ggplot2)

set.seed(1071)

# data generated as shown in the zoo vignette
Z.index <- as.Date(sample(12450:12500, 10))
Z.data <- matrix(rnorm(30), ncol = 3)
colnames(Z.data) <- c("Aa", "Bb", "Cc")
Z <- zoo(Z.data, Z.index)

tidy(Z)

ggplot(tidy(Z), aes(index, value, color = series)) +
  geom_line()
  
ggplot(tidy(Z), aes(index, value)) +
  geom_line() +
  facet_wrap(~ series, ncol = 1)

Zrolled <- rollmean(Z, 5)
ggplot(tidy(Zrolled), aes(index, value, color = series)) +
  geom_line()




cleanEx()
nameEx("tidy_optim")
### * tidy_optim

flush(stderr()); flush(stdout())

### Name: tidy_optim
### Title: Tidy a(n) optim object masquerading as list
### Aliases: tidy_optim optim_tidiers tidy.optim

### ** Examples


func <- function(x) {
    (x[1] - 2)^2 + (x[2] - 3)^2 + (x[3] - 8)^2
}

o <- optim(c(1, 1, 1), func)

tidy(o)
glance(o)




cleanEx()
nameEx("tidy_svd")
### * tidy_svd

flush(stderr()); flush(stdout())

### Name: tidy_svd
### Title: Tidy a(n) svd object masquerading as list
### Aliases: tidy_svd svd_tidiers

### ** Examples


mat <- scale(as.matrix(iris[, 1:4]))
s <- svd(mat)

tidy_u <- tidy(s, matrix = "u")
tidy_u

tidy_d <- tidy(s, matrix = "d")
tidy_d

tidy_v <- tidy(s, matrix = "v")
tidy_v

library(ggplot2)
library(dplyr)

ggplot(tidy_d, aes(PC, percent)) +
    geom_point() +
    ylab("% of variance explained")

tidy_u %>%
    mutate(Species = iris$Species[row]) %>%
    ggplot(aes(Species, value)) +
    geom_boxplot() +
    facet_wrap(~ PC, scale = "free_y")





cleanEx()
nameEx("tidy_xyz")
### * tidy_xyz

flush(stderr()); flush(stdout())

### Name: tidy_xyz
### Title: Tidy a(n) xyz object masquerading as list
### Aliases: tidy_xyz xyz_tidiers

### ** Examples


A <- list(x = 1:5, y = 1:3, z = matrix(runif(5 * 3), nrow = 5))
image(A)
tidy(A)




cleanEx()
nameEx("vector_tidiers")
### * vector_tidiers

flush(stderr()); flush(stdout())

### Name: tidy.numeric
### Title: Tidy atomic vectors
### Aliases: tidy.numeric tidy.character tidy.logical

### ** Examples


## Not run: 
##D x <- 1:5
##D names(x) <- letters[1:5]
##D tidy(x)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
