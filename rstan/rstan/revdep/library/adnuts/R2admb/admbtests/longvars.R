library(R2admb)
## reed frog example, but make variable names longer
## source("~/lib/R/pkgs/R2admb/pkg/R/admb-funs.R")

ReedfrogSizepred <- 
  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L))

setup_admb()
fn <- "longvars"
tplfile <- paste(fn,"tpl",sep=".")
file.copy(system.file("tplfiles",tplfile,package="R2admb"),
          tplfile)

m1 <- do_admb(fn,
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
              profile=TRUE,
              profile.opts=list(pars=c("ccc","dddd","ggggggg")),
              params=list(ccc=0.45,dddd=13,ggggggg=1),
              bounds=list(ccc=c(0,1),dddd=c(0,50),ggggggg=c(-1,25)),
              run.opts=run.control(checkparam="write",
                checkdata="write"),
              workdir="tmp",
              verbose=TRUE)



m1P <- do_admb(fn,
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
               params=as.list(coef(m1)),
##              bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
##                bounds=list(c=c(0.1,1),d=c(5,50),g=c(-1,25)),
               bounds=list(c=c(0.3,0.6),d=c(10,15),g=c(15,20)),
               run.opts=run.control(checkparam="write",
                 checkdata="write"),
               profile=TRUE,
               verbose=TRUE,
              profpars=c("c","d","g"))
## argh -- must fix this.  don't know where the errors are coming from --
##  was OK in a previous version

m1MC <- do_admb(fn,
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
                params=list(c=0.45,d=13,g=1),
                bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
                run.opts = run.control(checkparam="write",
                  checkdata="write"),
                mcmc=TRUE,
                mcmc.opts=mcmc.opts(mcmcpars=c("c","d","g")))
save("m1","m1P","m1MC",file="Reedfrog_runs.RData")

## clean up
unlink(c("reedfrogsizepred0","reedfrogsizepred0.tpl","reedfrogsizepred0_gen.tpl"))
