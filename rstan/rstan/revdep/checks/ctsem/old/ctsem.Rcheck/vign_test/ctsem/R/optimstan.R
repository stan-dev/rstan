#' Optimize / importance sample a stan or ctStan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init vector of unconstrained parameter values, or character string 'random' to initialise with 
#' random values very close to zero.
#' @param sampleinit either NA, or an niterations * nparams matrix of samples to initialise importance sampling.
#' @param deoptim Do first pass optimization using differential evolution? Slower, but better for cases with multiple 
#' minima / difficult optimization.
#' @param stochastic Logical. Use stochastic gradient descent instead of ucminf (bfgs with trust region) optimizer.
#' Generally more robust, a little slower with few parameters, faster with many.
#' @param plotsgd Logical. If TRUE, plot iteration details when using stochastic optimizer.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{\link[DEoptim]{DEoptim.control}}.
#' @param nopriors logical.f If TRUE, any priors are disabled -- sometimes desirable for optimization. 
#' @param tol absolute object tolerance
#' @param cores Number of cpu cores to use.
#' @param isloops Number of iterations of adaptive importance sampling to perform after optimization.
#' @param isloopsize Number of samples per iteration of importance sampling.
#' @param finishsamples Number of samples to use for final results of importance sampling.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent 
#' importance sampling, at risk of truncating tails.
#'
#' @return ctStanFit object
#' @importFrom ucminf ucminf
#' @importFrom Matrix bdiag
#' @importFrom utils head tail
#' @export
#'
#' @examples
#'  sunspots<-sunspot.year
#'  sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#'  id <- 1
#'  time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#' 
#' #setup model
#'  ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1,
#'   manifestNames='sunspots',
#'   latentNames=c('ss_level', 'ss_velocity'),
#'    LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#'    DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
#'    MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#'    MANIFESTVAR=diag(0,1),
#'    CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'    T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
#'    DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#' 
#'  ssmodel$pars$indvarying<-FALSE #Because single subject
#'  ssmodel$pars$offset[14]<- 44 #Because not mean centered
#'  ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality
#' 
#' #fit using optimization without importance sampling
#' ssfit <- ctStanFit(datalong[1:50,], #limited data for example
#'   ssmodel, optimize=TRUE,optimcontrol=list(deoptim=FALSE,isloops=0,finishsamples=50))
#' 
#' #output
#' summary(ssfit)
optimstan <- function(standata, sm, init='random',sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-12,
  decontrol=list(),
  stochastic = TRUE,
  plotsgd=FALSE,
  isloops=0, isloopsize=1000, finishsamples=500, tdf=50,
  verbose=0,nopriors=FALSE,cores=1){
  
  standata$verbose=as.integer(verbose)
  standata$nopriors=as.integer(nopriors)
  
  if(is.null(decontrol$steptol)) decontrol$steptol=5 
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  
  message('Optimizing...')
  
  betterfit<-TRUE
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    # if(nopriors){
    #   standata$nopriors <- 0
    #   suppressWarnings(suppressOutput(optimfit <- optimizing(sm,standata, hessian=FALSE, iter=400, init=0,as_vector=FALSE,
    #     tol_obj=1e-8, tol_rel_obj=0,init_alpha=.000001, tol_grad=0,tol_rel_grad=1e7,tol_param=1e-5,history_size=100),verbose=verbose))
    #   init=optimfit$par
    #   standata$nopriors <- 1
    # }
    
    # suppressMessages(suppressWarnings(suppressOutput(smf<-sampling(sm,iter=0,chains=0,init=0,data=standata,check_data=FALSE, 
    # control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
    smf <- stan_reinitsf(sm,standata)
    npars=get_num_upars(smf)
    if(all(init %in% 'random')) init <- rnorm(npars, 0, .001)
    if(all(init == 0)) init <- rep(0,npars)
    
    if(is.na(sampleinit[1])){
      
      if(deoptim){ #init with DE
        # require(DEoptim)
        if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
        
        decontrollist <- c(decontrol,DEoptim.control())
        decontrollist <- decontrollist[unique(names(decontrollist))]
        
        lp2 = function(parm) {
          out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
          if(class(out)=='try-error') {
            out=-1e200
          }
          return(-out)
        }
        
        deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
        deinit[2,] <- rnorm(npars,0,.0002)
        if(length(init)>1 & try2) {
          deinit[1,] <- unconstrain_pars(smf,init)
          if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
        }
        decontrollist$initialpop=deinit
        decontrollist$NP = NP
        optimfitde <- suppressWarnings(DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
          control = decontrollist))
        # init=constrain_pars(object = smf,optimfitde$optim$bestmem)
        init=optimfitde$optim$bestmem
      }
      
      
      # suppressWarnings(suppressOutput(optimfit <- optimizing(sm,standata, hessian=FALSE, iter=1e6, init=init,as_vector=FALSE,draws=0,constrained=FALSE,
      #   tol_obj=tol, tol_rel_obj=0,init_alpha=.001, tol_grad=0,tol_rel_grad=0,tol_param=0,history_size=50,verbose=verbose),verbose=verbose))
      
      gradout <- c()
      bestlp <- -Inf
      
      lp<-function(parm) {
        # print((parm))
        out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=-Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          if(out[1] > bestlp) {
            bestlp <<- out[1]
            gradout <<- attributes(out)$gradient
          }
        }
        return(-out[1])
      }
      
      grffromlp<-function(parm) {
        return(-gradout)
      }
      
      parbase=par()
      
      sgd <- function(init,stepbase=1e-4,nsubjects=1,gmeminit=0,gmemmax=.9,maxparchange = .1,
        minparchange=1e-16,maxiter=5000,perturbpercent=0,nconvergeiter=20, itertol=1e-3, deltatol=1e-5){
        pars=init
        bestpars = pars
        step=rep(stepbase,length(init))
        g=rnorm(length(init),0,.001)
        gsmooth=g
        gpersist=g
        signg=sign(g)
        oldsigng = g
        gmemory <- gmeminit
        oldgmemory <- gmemory
        oldlpdif <- 0
        lpdif <- 0
        maxlp <- -Inf
        i=0
        lp<-c()
        oldlp <- -Inf
        converged <- FALSE
        while(!converged && i < maxiter){
          i = i + 1
          accepted <- FALSE
          while(!accepted){
            whichpars = sample(1:length(pars), ceiling(length(pars)*perturbpercent))
            newpars = bestpars #* .5 + bestpars * .5
            parupdate =  exp(rnorm(length(step),0,step)) * step  * gsmooth #(sign(g) +  sign(g)*(abs(gsmooth) / mean(abs(gsmooth)))) #
            if(i %% 28 == 0) step[whichpars] = step[whichpars] * 10
            newpars = newpars + parupdate
            # subjects <- sample(unique(standata$subject),nsubjects)
            # standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,standata$ndatapoints) #
            # smf<-stan_reinitsf(sm,standata)
            lpg = try(log_prob(smf, upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
            if(class(lpg) !='try-error' && !is.nan(lpg[1])) accepted <- TRUE else step <- step * .1
            if(i < 20 && i > 1 && lpg[1] < lp[i-1]) {
              accepted <- FALSE
              step = step * .1
            }
          }
          lp=c(lp,lpg[1])
          pars <- newpars
          oldsigng=signg
          # gpersist <- gpersist * .99 + .01 * abs(attributes(lpg)$gradient -g) / (abs(attributes(lpg)$gradient)+abs(g)/2)
          # print(gpersist[1:10])
          g=attributes(lpg)$gradient
          oldgsmooth = gsmooth
          gsmooth= gsmooth*gmemory + (1-gmemory)*g
          signg=sign(gsmooth)
          # step=exp(mean(log(step))+(.99*(log(step)-mean(log(step)))))
          step[oldsigng == signg] = step[oldsigng == signg] * 1.1 #exp((1-gmemory)/2)
          step[oldsigng != signg] = step[oldsigng != signg] / 1.2 #exp((1-gmemory)/2)
          if(lp[i] >= maxlp) {
            step = step * 1.05 #exp((1-gmemory)/8)
            bestpars <- pars
            maxlp <- lp[i]
          } 
          if(i > 1 && lp[i] < lp[i-1]) {
            # signg <- oldsigng
            # gsmooth = oldgsmooth
            # pars <- bestpars
            step = step  / max( 1.5, (-10*(lp[i] - lp[i-1]) / sd(head(tail(lp,20),10)))) #exp((1-gmemory)/4)
          }
          # if(i %%10 ==0) gmemory = min(gmemory+.1,gmemmax)# * exp(mean(sign(diff(tail(lp,20)))))
          if(i %%20 ==0) gmemory =  max(gmeminit, min(gmemmax, 1.6*(1-(log(sd(tail(lp,20)) ) -log(itertol)) / (log(sd(head(lp,20)))-log(itertol)))* (1-gmeminit) + gmeminit))
          
          # if(i > 30 && i %% 10 == 0) {
          #   lpdif <- sum(diff(tail(lp,10)))
          #   oldlpdif <- sum(diff(head(tail(lp,10),20)))
          #   if(oldlpdif >= lpdif) gmemory <- oldgmemory
          #   proposal = gmemory*2-oldgmemory
          #   gmemory <- min(gmemmax, max(0, proposal + runif(1,-.2,.2)))
          #   oldgmemory <- gmemory
          #   if(i < 50 && sign(lpdif)==-1 & sign(oldlpdif)==-1) {
          #     message('diverging')
          #     # gmemory=gmemory*.8
          #   }
          #   oldgmemory <- gmemory
          # }
          
          step[step > maxparchange] <- maxparchange
          step[step < minparchange] <- minparchange
          
          if(plotsgd){
            par(mfrow=c(3,1))
            plot(pars)
            plot(log(step))
            plot(tail(log(abs(lp)),500),type='l')
            print(i)
            message(paste0('Iter = ',i, '   Best LP = ', maxlp,'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
          }
          
          #check convergence
          if(i > 30){
            if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
            if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
          }
        }
        return(list(itervalues = lp, value = maxlp,par=bestpars) )
      }
      # browser()
      # if(!deoptim & standata$nopriors == 0 ) init='random'
      
      if(stochastic=='auto' && npars > 50){
        message('> 50 parameters and stochastic="auto" so stochastic gradient descent used')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE
      
      if(!deoptim & standata$nopriors == 1 ){ #init using priors
        standata$nopriors <- as.integer(0)
        smf <- stan_reinitsf(sm,standata)
        if(!stochastic) optimfit <- ucminf(init,fn = lp,gr = grffromlp,control=list(grtol=1e-4,xtol=tol*1e4,maxeval=10000))
        if(stochastic) optimfit <- sgd(init, stepbase=1e-4,nsubjects=standata$nsubjects, gmemmax=ifelse(npars > 50, .95, .9)) 
        standata$nopriors <- as.integer(1)
        smf <- stan_reinitsf(sm,standata)
        init = optimfit$par #rstan::constrain_pars(object = smf, optimfit$par)
      }
      
      
      if(!stochastic) {
        optimfit <- ucminf(init,fn = lp,gr = grffromlp,control=list(grtol=1e-8,xtol=tol,maxeval=10000),hessian=2)
        init = optimfit$par
        ucminfcov <- optimfit$invhessian
      }
      # sink('fail.txt')
      optimfit <- sgd(init, nsubjects=standata$nsubjects,perturbpercent = .0) 
      # sink()
      # browser()
      # optimizing(object = sm, standata)
      est1=constrain_pars(smf,optimfit$par)
      bestfit <-optimfit$value
      # smf<-new(sm@mk_cppmodule(sm),standata,0L,rstan::grab_cxxfun(sm@dso))
      
      est2=optimfit$par #unconstrain_pars(smf, est1)
    }
    
    if(!estonly){
      lp<-function(parm) {
        out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
        if(class(out)=='try-error') {
          out=-Inf
        }
        return(out)
      }
      
      grf<-function(parm,...) {
        out=try(grad_log_prob(smf, upars=parm, adjust_transform = TRUE))
        if(class(out)=='try-error') {
          out=rep(NA,length(parm))
        }
        return(out)
      }
      
      grmat<-function(func,pars,step=1e-12){
        gradout<-matrix(NA,nrow=length(pars),ncol=length(pars))
        
        for(i in 1:length(pars)){
          stepsize <- step *10
          while((any(is.na(gradout[i,])) || gradout[i,i] >=0)  && stepsize > 1e-12){
            stepsize <- stepsize * .1
            uppars<-pars
            downpars<-pars
            uppars[i]<-pars[i]+stepsize
            downpars[i]<-pars[i]-stepsize
            gradout[i,]<-((func(uppars)) - (func(downpars)))/stepsize/2
          }
        }
        return(t(gradout))
      }
      
      # A more numerically stable way of calculating log( sum( exp( x ))) Source:
      # http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
      log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
      }
      
      
      if(is.na(sampleinit[1])){
          
          hess=grmat(func=grf,pars=est2,step=1e-8)
          if(any(is.na(hess))) stop(paste0('Hessian could not be computed for pars ', paste0(which(apply(hess,1,function(x) any(is.na(x))))), ' -- consider reparameterising.',collapse=''))
          hess = (hess/2) + t(hess/2)
          mchol=try(t(chol(solve(-hess))),silent=TRUE)
          if(class(mchol)=='try-error') {
            message('Hessian not positive-definite -- check importance sampling convergence with isdiag')
            npd <- TRUE
          } else npd <- FALSE
          # if(class(mchol)=='try-error') {
          mcov=MASS::ginv(-hess) #-optimfit$hessian)
          mcov=as.matrix(Matrix::nearPD(mcov)$mat)
      }
      
      if(!is.na(sampleinit[1])){
        mcov = cov(sampleinit)*1.5+diag(1e-6,ncol(sampleinit))
        est2 = apply(sampleinit,2,mean)
        bestfit = 9e100
        optimfit <- suppressWarnings(list(par=sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
      }
      
      mcovl <- list()
      mcovl[[1]]=mcov
      delta=list()
      delta[[1]]=est2
      samples <-c()
      resamples <- c()
      prop_dens <-c()
      target_dens<-c()
      sample_prob<-c()
      counter <- 0
      ess <- 0
      qdiag<-0
      
      cl <- parallel::makeCluster(cores, type = "PSOCK")
      parallel::clusterExport(cl, c('sm','standata'),environment())
      
      if(isloops == 0) {
        nresamples = finishsamples
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + t(chol(mcovl[[1]])) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        message('Importance sampling not done -- interval estimates via Hessian based sampling only')
      }
      
      if(isloops > 0){
        message('Adaptive importance sampling, loop:')
        j <- 0
        while(j < isloops){
          j<- j+1
          message(paste0('  ', j, ' / ', isloops, '...'))
          if(j==1){
            # if(!npd) 
            # browser()
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            # 
            # gensigstates <- function(t0means, t0chol,steps){
            #   out <- NA
            #   nlatent=nrow(t0chol)
            #   t0states <- matrix(t0means,byrow=TRUE,nlatent*2,nlatent)
            #   t0base <- matrix(t0means,byrow=TRUE,nlatent,nlatent)
            #   for(sqrtukfadjust in steps){
            #     sigpoints <- t(chol(as.matrix(Matrix::bdiag((t0chol%*%t(t0chol))))))*sqrtukfadjust
            #     t0states[1:(nlatent),] =  t0base + t(sigpoints)
            #     t0states[(1+nlatent):(nlatent*2),] = t0base - t(sigpoints)
            #     if(is.na(out[1])) out <- t0states else out <- rbind(out,t0states)
            #   }
            #   return(out)
            # }
            # samples=gensigstates(delta[[j]],mcovl[[j]],5000*(.1^(exp(seq(0,1,length.out=ceiling(isloopsize/npars/2))))))
            # samples=samples[sample(1:nrow(samples),isloopsize),]
            # 
            # if(npd){
            #   samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            #   prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize*10), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
            #   samples <- samples[prop_dens > quantile(prop_dens,probs = .9),]
            #   prop_dens <- prop_dens[prop_dens > quantile(prop_dens,probs = .9)]
            # }
          } else {
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = cov(resamples) #+diag(1e-12,ncol(samples))
            samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf))
          }
          # if(j > 1 || !npd) 
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          
          parallel::clusterExport(cl, c('samples'),environment())
          
          target_dens[[j]] <- unlist(parallel::parLapply(cl, parallel::clusterSplit(cl,1:isloopsize), function(x){
            eval(parse(text=paste0('library(rstan)')))
            
            smf <- stan_reinitsf(sm,standata)
            
            lp<-function(parm) {
              out<-try(log_prob(smf, upars=parm, adjust_transform = TRUE, gradient=FALSE),silent = TRUE)
              if(class(out)=='try-error') {
                out=-Inf
              }
              return(out)
            }
            out <- apply(tail(samples,isloopsize)[x,],1,lp)
            
            try(dyn.unload(file.path(tempdir(), paste0(smf@stanmodel@dso@dso_filename, .Platform$dynlib.ext))),silent = TRUE)
            return(out)
            
          }))
          
          if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
          if(any(target_dens[[j]] > bestfit && (j < isloops && !try2))){
            oldfit <- bestfit
            try2 <- TRUE
            bestfit<-max(target_dens[[j]],na.rm=TRUE)
            betterfit<-TRUE
            # init = rstan::constrain_pars(object = smf, samples[which(unlist(target_dens) == bestfit),])
            init = samples[which(unlist(target_dens) == bestfit),]
            message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
            break
          }
          nresamples = ifelse(j==isloops,finishsamples,5000)
          
          
          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # browser()
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)
          
          sample_prob <- c(sample_prob,exp((weighted_dens - log_sum_exp(weighted_dens)))) #sum to 1 for each iteration, normalise later
          # if(j==isloops) isloopsize = length(sample_prob) #on last loop use all samples for resampling
          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
          # points(target_dens2[sample_prob> (1/isloopsize * 10)], prop_dens[sample_prob> (1/isloopsize * 10)],col='red')
          resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE),
            prob = sample_prob / sum(sample_prob))
          # resample_i <- sample(tail(1:nrow(samples),isloopsize), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE), 
          #   prob = tail(sample_prob,isloopsize) / sum(tail(sample_prob,isloopsize) ))
          if(j < isloops){
            message(paste0(length(unique(resample_i)), ' unique samples drawn, from ', nresamples,' resamples of ', nrow(samples),' actual, probability sd = ', sd(sample_prob)))
            if(length(unique(resample_i)) < 100) {
              message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
              # return(est)
            }
          }
          resamples <- samples[resample_i, , drop = FALSE]
          # points(target_dens2[resample_i],prop_dens[resample_i],col='blue')
          # resamples=mcmc(resamples)
          
          ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
          qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
            (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) ) 
          })))
          
        }
      }
    }
  }#end while no better fit
  if(!estonly){
    if(isloops==0) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]
    
    # parallel::stopCluster(cl)
    message('Computing quantities...')
 
    # cl <- parallel::makeCluster(min(cores,chains), type = "PSOCK")
    parallel::clusterExport(cl, c('relistarrays','resamples','sm','standata','optimfit'),environment())
    
    # target_dens <- c(target_dens,
    transformedpars <- try(parallel::parLapply(cl, parallel::clusterSplit(cl,1:nresamples), function(x){
      require(ctsem)
      Sys.sleep(.1)
      smf <- stan_reinitsf(sm,standata)
      Sys.sleep(.1)
      # smf<-new(sm@mk_cppmodule(sm),standata,0L,rstan::grab_cxxfun(sm@dso))
      out <- list()
      skeleton=est1
      for(li in 1:length(x)){
        flesh = try(unlist(rstan::constrain_pars(smf, resamples[x[li],])))
        if(class(flesh) == 'try-error') {
          flesh <- unlist(skeleton)
          flesh[] <- NA
        }
        names(flesh) <- c()
        out[[li]] <-relistarrays(flesh, skeleton)
      }
      return(out)
    }))
    
    missingsamps <-sapply(1:length(transformedpars[[1]]), function(x) any(is.na(unlist(transformedpars[[1]][[x]]))))
    nasampscount <- sum(missingsamps) 
    
    transformedpars[[1]] <- transformedpars[[1]][!missingsamps]
    nresamples <- nresamples - nasampscount
    if(nasampscount > 0) message(paste0(nasampscount,' NAs generated during final sampling of ', finishsamples, '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
    
    transformedpars<-unlist(transformedpars,recursive = FALSE)
    
    
    
    tostanarray <- function(flesh, skeleton){
      skelnames <- names(skeleton)
      skelstruc <- lapply(skeleton,dim)
      count=1
      npars <- ncol(flesh)
      niter=nrow(flesh)
      out <- list()
      for(ni in skelnames){
        if(prod(skelstruc[[ni]])>0){
          if(!is.null(skelstruc[[ni]])){
            out[[ni]] <- array(flesh[,count:(count+prod(skelstruc[[ni]])-1)],dim = c(niter,skelstruc[[ni]]))
            count <- count + prod(skelstruc[[ni]])
          } else {
            out[[ni]] <- array(flesh[,count],dim = c(niter))
            count <- count + 1
          }
        }
      }
      return(out)
    }
    
    transformedpars=try(tostanarray(flesh = matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = est1))
    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if(class(sds)=='try-error') sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, lest)),
      unlist(constrain_pars(smf, est2)),
      unlist(constrain_pars(smf, uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    
    parallel::stopCluster(cl)
    
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2, rawposterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  suppressWarnings(do.call(par,parbase)) #reset par in case plots done
  return(stanfit)
}

