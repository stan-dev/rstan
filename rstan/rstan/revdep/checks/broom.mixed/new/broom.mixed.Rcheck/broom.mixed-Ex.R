pkgname <- "broom.mixed"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('broom.mixed')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("augment.ranef.mer")
### * augment.ranef.mer

flush(stderr()); flush(stdout())

### Name: augment.ranef.mer
### Title: Augmentation for random effects (for caterpillar plots etc.)
### Aliases: augment.ranef.mer

### ** Examples

if (require("lme4")) {
   load(system.file("extdata","lme4_example.rda",package="broom.mixed"))
   rr <- ranef(lmm1,condVar=TRUE)
   aa <- broom::augment(rr)
   ## Q-Q plot:
   if (require(ggplot2) && require(dplyr)) {
      g0 <- ggplot(aa,aes(estimate,qq,xmin=lb,xmax=ub))+
          geom_errorbarh(height=0)+
          geom_point()+facet_wrap(~variable,scale="free_x")
      ## regular caterpillar plot:
      g1 <- ggplot(aa,aes(estimate,level,xmin=lb,xmax=ub))+
         geom_errorbarh(height=0)+
         geom_vline(xintercept=0,lty=2)+
         geom_point()+facet_wrap(~variable,scale="free_x")
      ## emphasize extreme values
      aa2 <- group_by(aa,grp,level)
      aa3 <- mutate(aa2, keep=any(estimate/std.error>2))
      ## Update caterpillar plot with extreme levels highlighted
      ##  (highlight all groups with *either* extreme intercept *or*
      ##   extreme slope)
      ggplot(aa3, aes(estimate,level,xmin=lb,xmax=ub,colour=factor(keep)))+
         geom_errorbarh(height=0)+
         geom_vline(xintercept=0,lty=2)+
         geom_point()+facet_wrap(~variable,scale="free_x")+
         scale_colour_manual(values=c("black","red"), guide=FALSE)
   }
}



cleanEx()
nameEx("brms_tidiers")
### * brms_tidiers

flush(stderr()); flush(stdout())

### Name: brms_tidiers
### Title: Tidying methods for a brms model
### Aliases: brms_tidiers tidy.brmsfit glance.brmsfit augment.brmsfit

### ** Examples

 ## original model
 ## Not run: 
##D     brms_crossedRE <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
##D            iter = 500, chains = 2)
##D  
## End(Not run)
 if (require("brms")) {
   ## load stored object
   load(system.file("extdata", "brms_example.rda", package="broom.mixed"))

   fit <- brms_crossedRE
   tidy(fit)
   tidy(fit, parameters = "^sd_", conf.int = FALSE)
   tidy(fit, effects = "fixed", conf.method="HPDinterval")
   tidy(fit, effects = "ran_vals")
   tidy(fit, effects = "ran_pars", robust = TRUE)
   # glance method
   glance(fit)
   ## this example will give a warning that it should be run with
   ## reloo=TRUE; however, doing this will fail
   ## because the \code{fit} object has been stripped down to save space
   suppressWarnings(glance(fit, looic = TRUE, cores = 1))
   head(augment(fit))
}




cleanEx()
nameEx("fixef.MCMCglmm")
### * fixef.MCMCglmm

flush(stderr()); flush(stdout())

### Name: fixef.MCMCglmm
### Title: Extract fixed effects from an 'MCMCglmm' object
### Aliases: fixef.MCMCglmm

### ** Examples

## Not run: 
##D   # a simple MCMCglmm model
##D   data(PlodiaPO)
##D   m <- MCMCglmm(PO ~ 1, random= ~ FSfamily, data=PlodiaPO, verbose=FALSE)
##D 
##D   # only extract average fixed effects
##D   fixef(m, use = "mean")
##D 
##D   # histogram of posterior samples of fixed effects
##D   hist(fixef(m))
##D   # matches the mean
##D   rowMeans(fixef(m))
## End(Not run)



cleanEx()
nameEx("gamlss_tidiers")
### * gamlss_tidiers

flush(stderr()); flush(stdout())

### Name: gamlss_tidiers
### Title: Tidying methods for gamlss objects
### Aliases: gamlss_tidiers tidy.gamlss

### ** Examples

if (requireNamespace("gamlss", quietly = TRUE) &&
    requireNamespace("gamlss.data", quietly = TRUE)) {
    data(abdom, package="gamlss.data")
    ## Not run: 
##D          mod <- gamlss(y~pb(x), sigma.fo=~pb(x), family=BCT,
##D                        data=abdom, method=mixed(1,20))
##D     
## End(Not run)
    ## load stored object
    mod <- readRDS(system.file("extdata", "gamlss_example.rds",
                   package="broom.mixed"))
    tidy(mod)
}




cleanEx()
nameEx("glmmTMB_tidiers")
### * glmmTMB_tidiers

flush(stderr()); flush(stdout())

### Name: glmmTMB_tidiers
### Title: Tidying methods for glmmTMB models
### Aliases: glmmTMB_tidiers tidy.glmmTMB augment.glmmTMB glance.glmmTMB

### ** Examples

if (require("glmmTMB") && require("lme4")) {
    data("sleepstudy",package="lme4")
    ## original model:
    ## Not run: 
##D         lmm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
##D     
## End(Not run)
    ## load stored object
    load(system.file("extdata","glmmTMB_example.rda",package="broom.mixed"))
    tidy(lmm1)
    tidy(lmm1, effects = "fixed")
    tidy(lmm1, effects = "fixed", conf.int=TRUE)
    tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="uniroot")
    ## FIX: tidy(lmm1, effects = "ran_vals", conf.int=TRUE)
    head(augment(lmm1, sleepstudy))
    glance(lmm1)

    ## original model:
    ##  glmm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
    ##                  data = cbpp, family = binomial, weights=size)
    tidy(glmm1)
    tidy(glmm1, effects = "fixed")
    head(augment(glmm1, cbpp))
    head(augment(glmm1, cbpp, type.residuals="pearson"))
    glance(glmm1)
}



cleanEx()
nameEx("glmmadmb_tidiers")
### * glmmadmb_tidiers

flush(stderr()); flush(stdout())

### Name: glmmadmb_tidiers
### Title: Tidying methods for glmmADMB models
### Aliases: glmmadmb_tidiers glmmADMB_tidiers tidy.glmmadmb
###   augment.glmmadmb glance.glmmadmb

### ** Examples


if (require("glmmADMB") && require("lme4")) {
    ## original model
    ## Not run: 
##D         data("sleepstudy", package="lme4")
##D         lmm1 <- glmmadmb(Reaction ~ Days + (Days | Subject), sleepstudy,
##D                          family="gaussian")
##D     
## End(Not run)
    ## load stored object
    load(system.file("extdata","glmmADMB_example.rda",package="broom.mixed"))
    tidy(lmm1, effects = "fixed")
    tidy(lmm1, effects = "fixed", conf.int=TRUE)
    ## tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="profile")
    ## tidy(lmm1, effects = "ran_vals", conf.int=TRUE)
    head(augment(lmm1, sleepstudy))
    glance(lmm1)

    glmm1 <- glmmadmb(cbind(incidence, size - incidence) ~ period + (1 | herd),
                  data = cbpp, family = "binomial")
    tidy(glmm1)
    tidy(glmm1, effects = "fixed")
    head(augment(glmm1, cbpp))
    glance(glmm1)

}



cleanEx()
nameEx("lme4_tidiers")
### * lme4_tidiers

flush(stderr()); flush(stdout())

### Name: lme4_tidiers
### Title: Tidying methods for mixed effects models
### Aliases: lme4_tidiers tidy.merMod tidy.rlmerMod augment.merMod
###   glance.merMod

### ** Examples


if (require("lme4")) {
    ## original model
    ## Not run: 
##D         lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##D     
## End(Not run)
    ## load stored object
    load(system.file("extdata", "lme4_example.rda", package="broom.mixed"))
    tidy(lmm1)
    tidy(lmm1, effects = "fixed")
    tidy(lmm1, effects = "fixed", conf.int=TRUE)
    tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="profile")
    ## lmm1_prof <- profile(lmm1) # generated by extdata/runexamples
    tidy(lmm1, conf.int=TRUE, conf.method="profile", profile=lmm1_prof)
    ## conditional modes (group-level deviations from population-level estimate)
    tidy(lmm1, effects = "ran_vals", conf.int=TRUE)
    ## coefficients (group-level estimates)
    (rcoef1 <- tidy(lmm1, effects = "ran_coefs"))
    ## reconstitute standard coefficient-by-level table
    if (require(tidyr)) {
       spread(rcoef1,key=term,value=estimate)
    }
    head(augment(lmm1, sleepstudy))
    glance(lmm1)

    glmm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                  data = cbpp, family = binomial)
    tidy(glmm1)
    tidy(glmm1,exponentiate=TRUE)
    tidy(glmm1, effects = "fixed")
    ## suppress warning about influence.merMod
    head(suppressWarnings(augment(glmm1, cbpp)))
    glance(glmm1)

    startvec <- c(Asym = 200, xmid = 725, scal = 350)
    nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
                  Orange, start = startvec)
    ## suppress warnings about var-cov matrix ...
    op <- options(warn=-1)
    tidy(nm1)
    tidy(nm1, effects = "fixed")
    options(op)
    head(augment(nm1, Orange))
    glance(nm1)
    detach("package:lme4")
}
if (require("lmerTest")) {
   lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
   tidy(lmm1)
   glance(lmm1)
   detach("package:lmerTest")  # clean up
}



cleanEx()
nameEx("mcmc_tidiers")
### * mcmc_tidiers

flush(stderr()); flush(stdout())

### Name: tidy.MCMCglmm
### Title: Tidying methods for MCMC (Stan, JAGS, etc.) fits
### Aliases: tidy.MCMCglmm mcmc_tidiers tidyMCMC tidy.rjags tidy.stanfit
###   tidy.mcmc tidy.mcmc.list

### ** Examples

if (require("MCMCglmm")) {
  ## original model
  ## Not run: 
##D       mm0 <- MCMCglmm(Reaction ~ Days,
##D                  random = ~Subject, data = sleepstudy,
##D                  nitt=4000,
##D                  pr = TRUE
##D              )
##D    
## End(Not run)
   ## load stored object
   load(system.file("extdata","MCMCglmm_example.rda",
                                     package="broom.mixed"))
   tidy(mm0)
   tidy(mm1)
   tidy(mm2)
   tail(tidy(mm0,effects="ran_vals"))
}

# Using example from "RStan Getting Started"
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

model_file <- system.file("extdata", "8schools.stan", package = "broom.mixed")
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
## original model
## Not run: 
##D     set.seed(2015)
##D     rstan_example <- rstan::stan(file = model_file, data = schools_dat,
##D                          iter = 1000, chains = 2, save_dso = FALSE)
## End(Not run)
if (require(rstan)) {
   ## load stored object
   rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
   tidy(rstan_example)
   tidy(rstan_example, conf.int = TRUE, pars = "theta")
   td_mean <- tidy(rstan_example, conf.int = TRUE)
   td_median <- tidy(rstan_example, conf.int = TRUE, robust = TRUE)

  if (require(dplyr) && require(ggplot2)) {
    tds <- rbind(mutate(td_mean, method = "mean"),
             mutate(td_median, method = "median")) %>%
       mutate(type=ifelse(grepl("^theta",term),"theta",
            ifelse(grepl("^eta",term),"eta",
                  "other")))

     ggplot(tds, aes(estimate, term)) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),height=0) +
      geom_point(aes(color = method))+
      facet_wrap(~type,scale="free",ncol=1)
 } ## require(dplyr,ggplot2)
} ## require(rstan)
if (require(R2jags)) {
   ## see help("jags",package="R2jags")
   ## and  example("jags",package="R2jags")
   ## for details
   ## load stored object
   R2jags_example <- readRDS(system.file("extdata", "R2jags_example.rds", package = "broom.mixed"))
   tidy(R2jags_example)
   tidy(R2jags_example, conf.int=TRUE, conf.method="quantile")
}




cleanEx()
nameEx("nlme_tidiers")
### * nlme_tidiers

flush(stderr()); flush(stdout())

### Name: nlme_tidiers
### Title: Tidying methods for mixed effects models
### Aliases: nlme_tidiers tidy.lme augment.lme glance.lme tidy.gls

### ** Examples


if (require("nlme") && require("lme4")) {
    data("sleepstudy", package="lme4")
    ## original model
    ## Not run: 
##D          lmm1 <- lme(Reaction ~ Days, random=~ Days|Subject, sleepstudy)
##D     
## End(Not run)
    ## load stored object
    load(system.file("extdata","nlme_example.rda", package="broom.mixed"))
    tidy(lmm1)
    tidy(lmm1, effects = "fixed")
    tidy(lmm1, conf.int = TRUE)
    head(augment(lmm1, sleepstudy))
    glance(lmm1)

    startvec <- c(Asym = 200, xmid = 725, scal = 350)
    nm1 <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange,
                  fixed = Asym + xmid + scal ~1,
                  random = Asym ~1,
                  start = startvec)
    tidy(nm1)
    tidy(nm1, effects = "fixed")
    head(augment(nm1, Orange))
    glance(nm1)

    gls1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
                         correlation = corAR1(form = ~ 1 | Mare))
    tidy(gls1)
    glance(gls1)
}




cleanEx()
nameEx("paramNamesMCMCglmm")
### * paramNamesMCMCglmm

flush(stderr()); flush(stdout())

### Name: paramNamesMCMCglmm
### Title: Extract the parameter names from an 'MCMCglmm' object
### Aliases: paramNamesMCMCglmm
### Keywords: internal

### ** Examples

## Not run: 
##D   # a simple MCMCglmm model
##D   if (require(MCMCglmm)) {
##D      data(PlodiaPO)
##D      m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose=FALSE, pr=TRUE)
##D   }
##D   # extract the parameter names
##D   paramNamesMCMCglmm(m)
## End(Not run)



cleanEx()
nameEx("ranef.MCMCglmm")
### * ranef.MCMCglmm

flush(stderr()); flush(stdout())

### Name: ranef.MCMCglmm
### Title: Extract random effects from an 'MCMCglmm' object
### Aliases: ranef.MCMCglmm

### ** Examples

## Not run: 
##D   # a simple MCMCglmm model
##D   data(PlodiaPO)
##D   m <- MCMCglmm(PO ~ 1, random= ~ FSfamily, data=PlodiaPO, pr=TRUE, verbose=FALSE)
##D 
##D   # only extract average fixed effects
##D   head(ranef(m, use = "mean"))
##D 
##D   # histogram of posterior samples of fixed effects
##D   hist(ranef(m)[1, ])
##D   # matches the mean
##D   rowMeans(ranef(m)[1:6, ])
## End(Not run)



cleanEx()
nameEx("ranefLevels")
### * ranefLevels

flush(stderr()); flush(stdout())

### Name: ranefLevels
### Title: Extract the levels of factors used for random effects in
###   'MCMCglmm' objects
### Aliases: ranefLevels

### ** Examples

## Not run: 
##D   # a simple MCMCglmm model
##D   data(PlodiaPO)
##D   m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose=FALSE)
##D 
##D   # extract the random effects levels
##D   ranefLevels(m, PlodiaPO)
## End(Not run)



cleanEx()
nameEx("rstanarm_tidiers")
### * rstanarm_tidiers

flush(stderr()); flush(stdout())

### Name: rstanarm_tidiers
### Title: Tidying methods for an rstanarm model
### Aliases: rstanarm_tidiers tidy.stanreg glance.stanreg

### ** Examples


if (require("rstanarm")) {
## Not run: 
##D #'     ## original model
##D     fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
##D                       iter = 300, chains = 2)
##D   
## End(Not run)
## load example data
fit <- readRDS(system.file("extdata", "rstanarm_example.rds", package="broom.mixed"))

  # non-varying ("population") parameters
  tidy(fit, conf.int = TRUE, prob = 0.5)
  tidy(fit, conf.int = TRUE, conf.method = "HPDinterval", prob = 0.5)

  # hierarchical sd & correlation parameters
  tidy(fit, effects = "ran_pars")

  # group-specific deviations from "population" parameters
  tidy(fit, effects = "ran_vals")

  # glance method
   glance(fit)
  ## Not run: 
##D      glance(fit, looic = TRUE, cores = 1)
##D   
## End(Not run)
} ## if require("rstanarm")



cleanEx()
nameEx("stdranef")
### * stdranef

flush(stderr()); flush(stdout())

### Name: stdranef
### Title: Extract standard deviation of "random" effects from an
###   'MCMCglmm' object
### Aliases: stdranef

### ** Examples

## Not run: 
##D   # a simple MCMCglmm model
##D   data(PlodiaPO)
##D   PlodiaPO <- within(PlodiaPO, {
##D     PO2 <- cut(PO, quantile(PO, c(0, .33, .66, 1)))
##D     plate <- factor(plate)
##D   })
##D 
##D   m <- MCMCglmm(PO2 ~ 1, random = ~ FSfamily + plate,
##D     family = "ordinal", data = PlodiaPO,
##D     prior = list(
##D       R = list(V = 1, fix = 1),
##D       G = list(
##D         G1 = list(V = 1, nu = .002),
##D         G2 = list(V = 1, nu = .002)
##D       )
##D     ), verbose=FALSE, thin=1, pr=TRUE)
##D 
##D   # summary of the model
##D   summary(m)
##D 
##D   # examples of extracting standard deviations of
##D   # different random effects on the linear predictor metric
##D   # or after transformation to probabilities (only for ordinal)
##D   stdranef(m, which = list(1), type = "lp")
##D   stdranef(m, which = list(2), type = "lp")
##D   stdranef(m, which = list(1, 2, c(1, 2)), type = "lp")
##D   stdranef(m, type = "lp")
##D 
##D   ## error because no 3rd random effect
##D   #stdranef(m, which = list(1, 2, 3), type = "lp")
##D 
##D   stdranef(m, which = list("FSfamily", "plate"), type = "lp")
##D 
##D   # mean standard deviations on the probability metric
##D   # also the full distributions, if desired in the Data slot.
##D   res <- stdranef(m, type = "response")
##D   res$M # means
##D   hist(res$Data$FSfamily[, 1]) # histogram
## End(Not run)



cleanEx()
nameEx("tidy.TMB")
### * tidy.TMB

flush(stderr()); flush(stdout())

### Name: tidy.TMB
### Title: Tidying methods for TMB models
### Aliases: tidy.TMB

### ** Examples

if (require("TMB")) {
    runExample("simple",thisR=TRUE)
    class(obj) <- "TMB"
    tidy(obj,conf.int=TRUE,conf.method="wald")
    tidy(obj,conf.int=TRUE,conf.method="uniroot")
}



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
