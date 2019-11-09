## ----opts,echo=FALSE,message=FALSE---------------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE)

## ----do_admb_fake,eval=FALSE---------------------------------------------
#  fitted.model <- do_admb(fn,data,params,
#                          run.opts=run.control(checkparam="write",
#                          checkdata="write"))

## ----libs,message=FALSE--------------------------------------------------
library("R2admb")
library("ggplot2") ## for pictures
theme_set(theme_bw())  ## cosmetic
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))
library("bbmle")

## ----dat1----------------------------------------------------------------
ReedfrogSizepred <- 
  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L))

## ----mlefit--------------------------------------------------------------
m0 <- mle2(Kill~dbinom(c*((TBL/d)*exp(1-TBL/d))^g,size=10),
           start=list(c=0.45,d=13,g=1),data=ReedfrogSizepred,
           method="L-BFGS-B",
           lower=c(c=0.003,d=10,g=0),
           upper=c(c=0.8,d=20,g=60),
           control=list(parscale=c(c=0.5,d=10,g=1)))

## ----predvals------------------------------------------------------------
TBLvec = seq(9.5,36,length=100)
predfr <- 
  data.frame(TBL=TBLvec,
             Kill=predict(m0,newdata=data.frame(TBL=TBLvec)))

## ----fig1,echo=FALSE-----------------------------------------------------
g1  <- ggplot(ReedfrogSizepred,
              aes(x=TBL,y=Kill/10))+
    geom_point()+stat_sum(aes(size=..n..))+
    geom_smooth(method="loess")+
    labs(size="n",x="Size (total body length",
         y="Proportion killed")+
    coord_cartesian(ylim=c(-0.05,0.55))
startest <- stat_function(fun = function(x) { 0.45*((x/13)*exp(1-x/13)) },
                          lty=2,colour="red")
g1+startest+
    geom_line(data=predfr,colour="purple",lty=2)

## ----admbfit_getruns,echo=FALSE------------------------------------------
load_doc <- function(x) {
    ##    Define a bit of magic to load results from previous runs: 
    ## would like parent.frame(2) above instead of
    ## hacking around with global environment
    ## ... but doesn't work as I thought?
    load(system.file("doc",x,package="R2admb"),envir=.GlobalEnv) ##
}
zz <- load_doc("Reedfrog_runs.RData")

## ----setup_admb,eval=FALSE-----------------------------------------------
#  setup_admb()

## ----rfs_setup-----------------------------------------------------------
rfs_params <- list(c=0.45,d=13,g=1) ## starting parameters
rfs_bounds <- list(c=c(0,1),d=c(0,50),g=c(-1,25)) ## bounds
rfs_dat <- c(list(nobs=nrow(ReedfrogSizepred),
                  nexposed=rep(10,nrow(ReedfrogSizepred))),
             ReedfrogSizepred)

## ----admbfit_fake,eval=FALSE---------------------------------------------
#  m1 <- do_admb("ReedfrogSizepred0",
#                data=rfs_dat,
#                params=rfs_params,
#                bounds=rfs_bounds,
#                run.opts=run.control(checkparam="write",
#                checkdata="write",clean=FALSE))
#  unlink(c("reedfrogsizepred0.tpl",
#           "reedfrogsizepred0_gen.tpl",
#           "reedfrogsizepred0")) ## clean up leftovers

## ----rffit2,eval=FALSE---------------------------------------------------
#  rfs_params$g <- 2
#  m2 <- do_admb("ReedfrogSizepred1",
#                data=rfs_dat,
#                params=rfs_params,
#                bounds=rfs_bounds,
#                run.opts=run.control(checkparam="write",
#                checkdata="write"))

## ----basic---------------------------------------------------------------
m1

## ----coef----------------------------------------------------------------
coef(m1)

## ----summary-------------------------------------------------------------
summary(m1)

## ----vcov----------------------------------------------------------------
vcov(m1)

## ----others--------------------------------------------------------------
c(logLik(m1),deviance(m1),AIC(m1))

## ----profrun,eval=FALSE,tidy=FALSE---------------------------------------
#  m1P <- do_admb("ReedfrogSizepred0",
#                 data=c(list(nobs=nrow(ReedfrogSizepred),
#                 nexposed=rep(10,nrow(ReedfrogSizepred))),
#                 ReedfrogSizepred),
#                 params=rfs_params,
#                 bounds=rfs_bounds,
#                 run.opts=run.control(checkparam="write",
#                 checkdata="write"),
#                 profile=TRUE,
#                 workdir=".",
#                 profile.opts=list(pars=c("c","d","g")))

## ----mleprof,cache=TRUE--------------------------------------------------
m0prof <- profile(m0)

## ----profcalcs2,echo=FALSE,warning=FALSE---------------------------------
tmpf <- function(p,w="prof") {
  pp <- log(m1P$prof[[p]][[w]][,2])
  pp <- max(pp)-pp
  data.frame(param=p,z=sqrt(2*pp),
             par.vals.c=NA,par.vals.d=NA,par.vals.g=NA,
             focal=m1P$prof[[p]][[w]][,1])
}
quadf <- function(p) {
    m <- coef(m0)[p]
    se <- stdEr(m0)[p]
    pvec <- seq(m-3*se,m+3*se,length=31)
    data.frame(param=p,z=abs(pvec-m)/se,
                focal=pvec,
                par.vals.c=NA,par.vals.d=NA,par.vals.g=NA)
}
proflist <- do.call(rbind,lapply(list("c","d","g"),tmpf))
profnlist <- do.call(rbind,lapply(list("c","d","g"),tmpf,w="prof_norm"))
quadlist <- do.call(rbind,lapply(list("c","d","g"),quadf))
pdat <- rbind(cbind(as.data.frame(m0prof),method="mle2"),
              cbind(proflist,method="ADMB"),
              cbind(profnlist,method="ADMB_norm"),
              cbind(quadlist,method="Wald"))

## ----profpic,echo=FALSE,fig.width=8,fig.height=3.2,warning=FALSE---------
ggplot(pdat,aes(x=focal,y=abs(z),group=method,colour=method))+geom_line()+
      geom_point(alpha=0.5)+
  facet_grid(.~param,scale="free_x")+ylim(0,3)+xlab("")+
      ylab(expression(Delta(sqrt(-2*L))))+
      geom_hline(yintercept=1.96,lty=2)+zmargin

## ----admbfakemc,eval=FALSE,tidy=FALSE------------------------------------
#  m1MC <- do_admb("ReedfrogSizepred0",
#                  data=rfs_dat,
#                  params=rfs_params,
#                  bounds=rfs_bounds,
#                  run.opts=run.control(checkparam="write",
#                    checkdata="write"),
#                  mcmc=TRUE,
#                  mcmc.opts=mcmc.control(mcmcpars=c("c","d","g")))
#  ## clean up leftovers:
#  unlink(c("reedfrogsizepred0.tpl",
#           "reedfrogsizepred0_gen.tpl",
#           "reedfrogsizepred0"))

## ----mchistplot,message=FALSE--------------------------------------------
plot(m1MC$hist)

## ----coda----------------------------------------------------------------
library("coda")
mmc <- as.mcmc(m1MC$mcmc)

## ----mctraceplot---------------------------------------------------------
library("lattice")
xyplot(mmc)

## ----rgdiags-------------------------------------------------------------
raftery.diag(mmc)
geweke.diag(mmc)

## ----gdiag,echo=FALSE----------------------------------------------------
gd <- geweke.diag(mmc)
gd1 <- gd[["z"]][1]

## ----effsize-------------------------------------------------------------
effectiveSize(mmc)

## ----hpd-----------------------------------------------------------------
HPDinterval(mmc)

## ----mcdensplot----------------------------------------------------------
densityplot(mmc)

## ----lme4,message=FALSE--------------------------------------------------
library(lme4)
if (as.numeric(R.version$major)<3) {
## FIXME
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 family = binomial, data = cbpp)
}

## ----toymats-------------------------------------------------------------
X <- model.matrix(~period,data=cbpp)
Zherd <- model.matrix(~herd-1,data=cbpp)

## ----toydat--------------------------------------------------------------
tmpdat <- list(X=X,Zherd=Zherd,
                 incidence=cbpp$incidence,size=cbpp$size,
                 nobs=nrow(cbpp))

## ----loadtoy1,echo=FALSE-------------------------------------------------
zz2 <- load_doc("toy1_runs.RData")

## ----fakerun2,eval=FALSE-------------------------------------------------
#  d1 <- do_admb("toy1",
#                data=tmpdat,
#                params=list(beta=rep(0,ncol(X)),sigma_herd=0.1),
#                bounds=list(sigma_herd=c(0.0001,20)),
#                re=list(u_herd=ncol(Zherd)),
#                run.opts=run.control(checkdata="write",checkparam="write"),
#                mcmc=TRUE,
#                mcmc.opts=mcmc.control(mcmc=20,mcmcpars=c("beta","sigma_herd")))

## ----testprofinput,echo=FALSE,eval=FALSE---------------------------------
#  do_admb("toy1", data=tmpdat,
#          params=list(beta=rep(0,ncol(X)),sigma_herd=0.1),
#          bounds=list(sigma_herd=c(0.0001,20)),
#          re=list(u_herd=ncol(Zherd)),
#          run.opts=run.control(checkdata="write",checkparam="write",
#            clean=FALSE))
#  run_admb("toy1_gen",profile=TRUE)
#  read_admb("toy1_gen",profile=TRUE)

## ----coefsumlmer---------------------------------------------------------
## FIXME
## coef(summary(gm1))

## ----coefsumadmb---------------------------------------------------------
coef(summary(d1))[1:5,]

## ----ranefs--------------------------------------------------------------
## FIXME
## plot(ranef(gm1)$herd[,1],coef(d1)[6:20]*coef(d1)["sigma_herd"],
##     xlab="glmer estimate",ylab="ADMB estimate")
## abline(a=0,b=1)

## ----mcmcconf------------------------------------------------------------
detach("package:lme4") ## HPDinterval definition gets in the way
HPDinterval(as.mcmc(d1$mcmc[,6:20]))

