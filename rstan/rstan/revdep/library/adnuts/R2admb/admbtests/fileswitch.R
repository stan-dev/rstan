library(R2admb)
## source("~/R/pkgs/r2admb/pkg/R/do_admb.R") ## newest version
setup_admb()

## setup
library(lme4) ## for cbpp data set
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                family = binomial, data = cbpp) 
X <- model.matrix(~period,data=cbpp)
Zherd <- model.matrix(~herd-1,data=cbpp)

tmpdat <- list(X=X,Zherd=Zherd,
                 incidence=cbpp$incidence,size=cbpp$size,
                 nobs=nrow(cbpp))

testfun <- function(fn) {
  do_admb(fn,
          data=tmpdat,
          re=TRUE,
          params=list(beta=rep(0,ncol(X)),sigma_herd=0.1),
          bounds=list(sigma_herd=c(0.0001,20)),
          re_vectors=c(u_herd=ncol(Zherd)),
          checkdata="write",checkparam="write")
}

file.copy("toy1.tpl","toy1_orig.tpl")
t1 <- testfun("toy1")
file.copy("toy1.tpl","TOY1.tpl")
t2 <- testfun("TOY1")
unlink("toy1.tpl")
t2 <- testfun("TOY1")

## clean up
file.copy("toy1_orig.tpl","toy1.tpl",overwrite=TRUE)
unlink("toy1_orig.tpl")
