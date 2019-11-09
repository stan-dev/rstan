
ctFitAuto <- function(datalong,manifestNames,idcol='id',timecol='time',
  sequence=c('addtrend','extendtrend','addbetween','extendbetween','adddynamics','extenddynamics','removetrend') ){
  
  trendwbw<-FALSE #trend includes bw difs at same time
  
  Tpoints=max(unlist(lapply(unique(datalong[,idcol]),function(x) 
    length(datalong[datalong[,idcol]==x, idcol]) )))
  
  out=list()
  
  if(!all(sequence %in% c('addtrend','extendtrend','addbetween',
    'extendbetween','adddynamics','extenddynamics','removetrend','removedynamics'))) stop('Invalid sequence!')
  
  for(mani in 1:length(manifestNames)){
    message(paste0('Fitting model to variable: ',manifestNames[mani]))
    out[[manifestNames[mani]]] <- list()
    modelcount = 0
    finished = FALSE
    
    #set empty current model
    cm<-list(n.latent=1)
    
    sequenceposition=0
    
    while(sequenceposition <=length(sequence)){
      modelcount = modelcount + 1
      out[[manifestNames[mani]]][[modelcount]] <- list()
      
      #reset next model addition
      nm<-list(addtrend=FALSE,addbetween=FALSE,adddynamics=FALSE,
        extendtrend=FALSE,extendbetween=FALSE,extenddynamics=FALSE,removetrend=FALSE,
        removedynamics=FALSE)

      if(modelcount > 1) {
        cm = bm #set current model to best model
        nm[[sequence[sequenceposition] ]] <- TRUE
      }
      
      
      if(trendwbw & (nm$addtrend | nm$extendtrend)) nm$extendbetween<-TRUE #disable for distinct trend and bw...
      

      if(modelcount==1){#configure base model
        message('Fitting with measurement error only')
        LAMBDA=matrix(1)
        DRIFT=matrix(-10,cm$n.latent,cm$n.latent)
        DIFFUSION=diag(.00001,cm$n.latent,cm$n.latent)
        T0MEANS=matrix(0,1,cm$n.latent)
        T0VAR=matrix(.00001)
        
        cm$ndynamics=0
        cm$dynlatents=c()
        cm$ntrend=0
        cm$trendlatents = c()
        cm$nbetween <- 0
        cm$betweenlatents=c()
      }
      
      if(nm$adddynamics){
        message('Adding stochastic process')
        cm$ndynamics=cm$ndynamics+1
        
        
        LAMBDA=bm$m$LAMBDA
        DRIFT=bm$m$DRIFT
        DIFFUSION=bm$m$DIFFUSION
        T0MEANS=bm$m$T0MEANS
        T0VAR=bm$m$T0VAR
        
        if(cm$ndynamics==1){
          cm$dynlatents=1
        DRIFT[1,1]=paste0('drift_dyn_',cm$dynlatents,cm$dynlatents)
        DIFFUSION[1,1]=paste0('diffusion_dyn_',cm$dynlatents,cm$dynlatents)
        T0VAR[1,1]=0
        }
        
        if(cm$ndynamics>1){
          cm$n.latent=cm$n.latent+1
          cm$dynlatents=c(cm$dynlatents,cm$n.latent)
          LAMBDA=cbind(bm$m$LAMBDA,1) 
          
          DRIFT=cbind( rbind(bm$m$DRIFT, 0), 0)
          DRIFT[cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]] <- paste0('drift_dyn',cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics])
          
          DIFFUSION=cbind(  rbind(cm$m$DIFFUSION,0) , c(rep(0,nrow(cm$m$DIFFUSION)),paste0('diffusion_dyn_',cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]))) #should this be on 1st or last order state?
          T0MEANS=rbind(bm$m$T0MEANS,0)
          
          T0VAR=cbind( rbind(bm$m$T0VAR, 0), 0) #stationary labels added after model spec
        }
        
      }
      
      if(nm$addbetween){ #add between subject differences
        message('Adding stable between subject differences')
        cm$n.latent=cm$n.latent+1
        cm$nbetween <- 1
        cm$betweenlatents=cm$n.latent
        LAMBDA=cbind(cm$m$LAMBDA,1) 
        DRIFT=cbind(  rbind(cm$m$DRIFT,0) , c(rep(0,cm$n.latent-1),-.00001))
        DIFFUSION=cbind(  rbind(cm$m$DIFFUSION,0) , c(rep(0,cm$n.latent-1),.00001))
        T0MEANS=rbind(cm$m$T0MEANS,0)
        T0VAR=cbind( rbind(bm$m$T0VAR,0) , c(rep(0,cm$n.latent-1),'traitvar')) 
      }
      
      if(nm$addtrend) {
        message('Adding deterministic trend')
        cm$n.latent=cm$n.latent+1
        cm$ntrend=1
        cm$trendlatents=cm$n.latent
        LAMBDA=cbind(bm$m$LAMBDA,1) 
        DRIFT=cbind(  rbind(bm$m$DRIFT,0) , c(rep(0,nrow(bm$m$DRIFT)),paste0('drift_trend',cm$trendlatents[cm$ntrend],cm$trendlatents[cm$ntrend])))
        DIFFUSION=cbind(  rbind(bm$m$DIFFUSION,0) , c(rep(0,nrow(bm$m$DIFFUSION)),.00001))
        T0MEANS=rbind(bm$m$T0MEANS,paste0('t0_trend_',cm$trendlatents))
        T0VAR=cbind( rbind(bm$m$T0VAR,0) , c(rep(0,nrow(bm$m$T0VAR)),0) )
      }
      
      if(nm$extendtrend) {
        message('Increasing order of deterministic trend')
        cm$n.latent=cm$n.latent+1
        cm$ntrend=cm$ntrend+1
        cm$trendlatents=c(cm$trendlatents,cm$n.latent)
        LAMBDA=cbind(bm$m$LAMBDA,0) 
        T0VAR=cbind( rbind(bm$m$T0VAR, 0), 0)
        
        DRIFT=cbind( rbind(bm$m$DRIFT, 0), 0)
        DRIFT[cm$trendlatents[cm$ntrend], cm$trendlatents[cm$ntrend-1]] <- paste0('drift_trend',cm$trendlatents[cm$ntrend],cm$trendlatents[cm$ntrend-1])
        DRIFT[cm$trendlatents[cm$ntrend-1],cm$trendlatents[cm$ntrend]] <- 1
        DRIFT[cm$trendlatents[cm$ntrend],cm$trendlatents[cm$ntrend]] <-paste0('drift_trend',cm$trendlatents[cm$ntrend],cm$trendlatents[cm$ntrend])
        DRIFT[cm$trendlatents[cm$ntrend-1],cm$trendlatents[cm$ntrend-1]]<- -.0001
      
        DIFFUSION=cbind(  rbind(bm$m$DIFFUSION,0) , c(rep(0,nrow(bm$m$DIFFUSION)),.0001))
        T0MEANS=rbind(bm$m$T0MEANS,paste0('t0_trend_',cm$ntrend))
      }
      
      
      if(nm$removetrend){
        message('Decreasing order of deterministic trend')
        LAMBDA=bm$m$LAMBDA[,- cm$trendlatents[cm$ntrend],drop=FALSE ]
        T0VAR=bm$m$T0VAR[- cm$trendlatents[cm$ntrend] ,- cm$trendlatents[cm$ntrend] ,drop=FALSE]
        
        DRIFT=bm$m$DRIFT[- cm$trendlatents[cm$ntrend] ,- cm$trendlatents[cm$ntrend] ,drop=FALSE]
       
        DIFFUSION=bm$m$DIFFUSION[- cm$trendlatents[cm$ntrend] ,- cm$trendlatents[cm$ntrend] ,drop=FALSE]
        T0MEANS=bm$m$T0MEANS[- cm$trendlatents[cm$ntrend] , ,drop=FALSE]
        
        cm$n.latent=cm$n.latent-1
        cm$trendlatents=cm$trendlatents[-cm$ntrend]
        cm$ntrend=cm$ntrend-1
        cm$dynlatents[-1] <- cm$dynlatents[-1] -1 #accounting for one less latent, all dynamics happen after trends (except first)
      }
      
      if(nm$removedynamics){
        message('Decreasing order of stochastic process')

        LAMBDA=bm$m$LAMBDA[,- cm$dynlatents[cm$ndynamics],drop=FALSE ]
        T0VAR=bm$m$T0VAR[- cm$dynlatents[cm$ndynamics] ,- cm$dynlatents[cm$ndynamics] ,drop=FALSE]
        
        DRIFT=bm$m$DRIFT[- cm$dynlatents[cm$ndynamics] ,- cm$dynlatents[cm$ndynamics] ,drop=FALSE]
        
        DIFFUSION=bm$m$DIFFUSION[- cm$dynlatents[cm$ndynamics] ,- cm$dynlatents[cm$ndynamics] ,drop=FALSE]
        T0MEANS=bm$m$T0MEANS[- cm$dynlatents[cm$ndynamics] , ,drop=FALSE]
        
        cm$trendlatents[cm$trendlatents > cm$dynlatents[cm$ndynamics] ] <- 
          cm$trendlatents[cm$trendlatents > cm$dynlatents[cm$ndynamics] ] -1 #accounting for one less latent
        cm$n.latent=cm$n.latent-1
        cm$dynlatents=cm$dynlatents[-cm$ndynamics]
        cm$ntrend=cm$ntrend-1

      }
      
      if(nm$extendbetween){
        message('Adding further between subjects differences')
        cm$nbetween <- cm$nbetween +1
        cm$betweenlatents <- c(cm$betweenlatents,cm$trendlatents[cm$nbetween-1]) #-1 because not adding  new latents but referencing trend
        
        if(!trendwbw){
        LAMBDA=bm$m$LAMBDA
        T0VAR=bm$m$T0VAR
        DRIFT=bm$m$DRIFT
        DIFFUSION=bm$m$DIFFUSION
        T0MEANS=bm$m$T0MEANS
        }
        
        # for(i in 1:cm$nbetween){ #covariance between all trend latents
        #   T0VAR[cm$betweenlatents[i:cm$nbetween], cm$betweenlatents[i]] <- paste0('t0var_between_',cm$betweenlatents[i:cm$nbetween],cm$betweenlatents[i])
        # }
        
        #covariance between
        for(i in cm$betweenlatents){
          for(j in cm$betweenlatents){
            if(i >=j) T0VAR[i,j] <- paste0('t0var_between_',i,j)
          }
        }

      }
      
      if(nm$extenddynamics) {
        message('Increasing order of stochastic process')
        cm$n.latent=cm$n.latent+1
        cm$ndynamics=cm$ndynamics+1
        cm$dynlatents=c(cm$dynlatents,cm$n.latent)
        LAMBDA=cbind(cm$m$LAMBDA,0) 
        
        DRIFT=cbind( rbind(bm$m$DRIFT, 0), 0)
        DRIFT[cm$dynlatents[cm$ndynamics], cm$dynlatents[cm$ndynamics-1]] <- paste0('drift_dyn',cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics-1])
        DRIFT[cm$dynlatents[cm$ndynamics-1],cm$dynlatents[cm$ndynamics]] <- 1
        DRIFT[cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]] <- paste0('drift_dyn',cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics])
        DRIFT[cm$dynlatents[cm$ndynamics-1],cm$dynlatents[cm$ndynamics-1]] <- -.00001
        
        # DIFFUSION=cbind(  rbind(cm$m$DIFFUSION,0) , c(rep(0,nrow(cm$m$DIFFUSION)),.00001)) #should this be on 1st or last order state?
        DIFFUSION=cbind(  rbind(cm$m$DIFFUSION,0) , c(rep(0,nrow(cm$m$DIFFUSION)),0))
        DIFFUSION[cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]] <-paste0('diffusion_dyn_',cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics])
        DIFFUSION[cm$dynlatents[cm$ndynamics-1],cm$dynlatents[cm$ndynamics-1]] <- 0
        
        T0MEANS=rbind(bm$m$T0MEANS,0)
        
        T0VAR=cbind( rbind(bm$m$T0VAR, 0), 0) #stationary labels added after model spec
      }
      
      
      #build ct model
      
      cm$m <- ctModel(n.manifest=1,n.latent=cm$n.latent,Tpoints=Tpoints, 
        manifestNames=manifestNames[mani],
        LAMBDA=LAMBDA,
        DRIFT=DRIFT,
        DIFFUSION=DIFFUSION,
        T0VAR=T0VAR,
        T0MEANS=T0MEANS)
      
      #generate omx model and do omx based changes
      
      if(modelcount == 1) fi=ctFit(dat = datalong,ctmodelobj = cm$m,dataform = 'long',
        objective = 'Kalman',meanIntervals = TRUE,verbose=0,fit=FALSE)
      
      
      if(modelcount > 1)  { #set inits to match lower order
        olddrift <- out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj$DRIFT$values
        oldT0var <- out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj$T0VARbase$values 
        olddiffusion <- out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj$DIFFUSIONbase$values 
        
        fi=ctFit(dat = datalong,ctmodelobj = cm$m,dataform = 'long',
        omxStartValues = omxGetParameters(out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj),
        carefulFit = FALSE,
        objective = 'Kalman',meanIntervals = TRUE,verbose=0,retryattempts = 1,fit=FALSE)
        
      }
      
      if(cm$ndynamics > 0){
     #set t0var for dynamics to stationary requires fixed parameters
        fi$mxobj <- mxModel(fi$mxobj,'T0VAR',remove = TRUE)
        fi$mxobj <- mxModel(fi$mxobj,
          mxAlgebra(T0VARchol %*% t(T0VARchol),name='T0VARalg'),
          mxMatrix(name='T0VAR',free=FALSE,values=0,nrow=cm$n.latent,ncol=cm$n.latent,
            labels=paste0('T0VARalg[',1:cm$n.latent,',',rep(1:cm$n.latent,each=cm$n.latent),']'))
        )
        
        for(i in 1:cm$ndynamics){ #covariance between all dynamic latents
          fi$mxobj$T0VAR$free[cm$dynlatents[i:cm$ndynamics], cm$dynlatents[i]] <- FALSE
          fi$mxobj$T0VAR$labels[cm$dynlatents[i:cm$ndynamics], cm$dynlatents[i]] <- paste0('asymDIFFUSION[',cm$dynlatents[i:cm$ndynamics],',',cm$dynlatents[i],']')
        }
      }
      
      
        if(nm$adddynamics){
          if(cm$ndynamics==1) fi$mxobj$DRIFT$values[1,1] <- olddrift[1,1]
          if(cm$ndynamics>1) fi$mxobj$DRIFT$values[cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]] <- -5
        }
        
        if(nm$extendtrend){
          #DRIFT
          fi$mxobj$DRIFT$values[cm$n.latent,cm$trendlatents[cm$ntrend-1]] <- 10 * 
            olddrift[cm$trendlatents[cm$ntrend-1],cm$trendlatents[cm$ntrend-1]]
          fi$mxobj$DRIFT$values[cm$n.latent,cm$n.latent] <- -10
          
          # #T0VAR
          # fi$mxobj$T0VARbase$values[cm$trendlatents[cm$ntrend],cm$trendlatents[1:(cm$ntrend-1)]] <- c(
          #   oldT0var[cm$trendlatents[cm$ntrend-1],cm$trendlatents[min(1,(cm$ntrend-2)):(cm$ntrend-2)]] , 1)
          # 
          # fi$mxobj$T0VARbase$values[cm$trendlatents[cm$ntrend],cm$trendlatents[cm$ntrend]] <- 0
        }

        
        if(nm$extenddynamics){  
          fi$mxobj$DRIFT$values[cm$n.latent,cm$dynlatents[cm$ndynamics-1]] <- 10 * 
            olddrift[cm$dynlatents[cm$ndynamics-1],cm$dynlatents[cm$ndynamics-1]]
          fi$mxobj$DRIFT$values[cm$n.latent,cm$n.latent] <- -10
          
          # fi$mxobj$DIFFUSIONbase$values[cm$dynlatents[cm$ndynamics],cm$dynlatents[cm$ndynamics]] <-  
          #   olddiffusion[cm$dynlatents[cm$ndynamics-1],cm$dynlatents[cm$ndynamics-1]]
          
          fi$mxobj$T0VARbase$values[cm$dynlatents[1],cm$dynlatents[cm$ndynamics-1]] <- ifelse(
            cm$ndynamics > 2, 
            oldT0var[cm$dynlatents[1],cm$dynlatents[cm$ndynamics-1]],
            1)
          fi$mxobj$T0VARbase$values[cm$n.latent,cm$n.latent] <- 0
        }
        
        #parameter restrictions
        indices <- matrix(NA,nrow=0,ncol=2)
        for(i in c(cm$trendlatents,cm$dynlatents)){
          for(j in c(cm$trendlatents,cm$dynlatents)){
            if( j <= i) indices <- rbind(indices,c(i,j))
          }
        }
        fi$mxobj$DRIFT$ubound[indices] <- 0
        
        
        fi$mxobj=mxRun(fi$mxobj)
      
      
      # #check misspecification via errors
      # kal<-list()
      # for(idi in unique(datalong[,idcol])){
      # kal[[idi]]=Kalman(kpars = ctModelFromFit(fi),
      #   datalong = datalong[datalong[,idcol]==idi, c(idcol,timecol,manifestNames[mani])],
      #   manifestNames = manifestNames[mani],
      #   latentNames=paste0('eta',1:n.latent))
      # }
      # 
      # # plot(kal[[1]]$time, kal[[1]]$y,pch=2)
      # # points(kal[[1]]$time, kal[[1]]$prederror,pch=2,type='l')
      # 
      # #error data file
      # errdat <- cbind(datalong[,idcol],
      #   datalong[,timecol],
      #   cbind(unlist(lapply(unique(datalong[,idcol]),function(x)
      #     kal[[x]]$yupd-kal[[x]]$y)))
      # )
      # colnames(errdat) = c('id','time','prederror')
      # 
      # errdat[,'id']=1 #remove this to test within subject model only?
      # errdat[,'time']=1:(nrow(errdat)) #and this, to account for time (but only makes sense when above is removed)
      # 
      # errmodelbase <- ctModel(n.manifest=1,n.latent=1,LAMBDA=diag(1),MANIFESTVAR=diag(0.00001,1),
      #   Tpoints=nrow(errdat),manifestNames = 'prederror',
      #   MANIFESTMEANS=diag(0,1),
      #   DRIFT=matrix('dr11'))
      # 
      # errmodeltest <- errmodelbase
      # errmodeltest$DRIFT[1,1]= -10
      # 
      # errmodelbasefit <-ctFit(dat=errdat,ctmodelobj = errmodelbase,stationary = c('T0VAR','T0MEANS'),
      #   objective = 'Kalman',dataform = 'long',meanIntervals = TRUE,verbose=2,retryattempts = 0)
      # 
      # errmodeltestfit <-ctFit(dat=errdat,ctmodelobj = errmodeltest,stationary = c('T0VAR','T0MEANS'),
      #   objective = 'Kalman',dataform = 'long',meanIntervals = TRUE,verbose=2,retryattempts = 0)
      
      
      
      #prepare output
      out[[manifestNames[mani]]][[modelcount]]$fit <- fi
      # out[[manifestNames[mani]]][[modelcount]]$misspec <- mxCompare(errmodelbasefit$mxobj,errmodeltestfit$mxobj)
      
      
      #compare model performance and prepare for next model (or finish)
      if(modelcount==1) {
        sequenceposition <- 1 #ready to start sequence now
        bestmodelindex <- 1
        bm <- cm
      }
      if(modelcount > 1){
        # compared<-mxCompare(out[[manifestNames[mani]]][[modelcount]]$fit$mxobj,
        #   out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj)
        
        aic <- 2 * length(out[[manifestNames[mani]]][[modelcount]]$fit$mxobj$estimate) + 
          out[[manifestNames[mani]]][[modelcount]]$fit$mxobj$output$fit
        
        aicbest <- 2 * length(out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj$estimate) + 
          out[[manifestNames[mani]]][[bestmodelindex]]$fit$mxobj$output$fit
        
        print(sequence[sequenceposition])
        print(paste0('AIC improvement = ', aicbest - aic))
        
        # if(compared$diffLL[2] < .1) browser()
        
        if(aic < aicbest -2) { #if model improved
          bm <- cm
          bestmodelindex=modelcount
          
          if(nm$addtrend | nm$adddynamics | nm$addbetween) { #since these are not recursive, go to next step even if successful
            sequenceposition=sequenceposition+1
          }
          
        }
        
        if(aic >= aicbest -2) {#if model didn't improve
          
          sequenceposition <- sequenceposition + 1
          if(nm$adddynamics & sequence[sequenceposition]=='extenddynamics') sequenceposition <- sequenceposition + 1 #if added dynamics was no good, can't extend
          
        }
      }
      if(sequenceposition <= length(sequence)){ #then check if the next step is viable
      if(sequence[sequenceposition]=='extendbetween' & bm$nbetween > bm$ntrend) sequenceposition <- sequenceposition +1 #if no trend to add differences for then skip
      if(sequence[sequenceposition]=='extenddynamics' & bm$ndynamics==0) sequenceposition <- sequenceposition +1 #if no dynamics to extend then skip
      if(sequence[sequenceposition]=='extendtrend' & bm$ntrend==0) sequenceposition <- sequenceposition +1 #if no trend to extend then skip
      if(sequence[sequenceposition]=='removetrend' & bm$ntrend==0) sequenceposition <- sequenceposition +1 #if no trend to remove then skip
      if(sequence[sequenceposition]=='removedynamics' & bm$ndynamics==0) sequenceposition <- sequenceposition +1 #if no dynamics to remove then skip
      }
      
      
    }
    out[[manifestNames[mani]]]$best <- out[[manifestNames[mani]]][[bestmodelindex]]$fit
  }
  return(out)
}
















