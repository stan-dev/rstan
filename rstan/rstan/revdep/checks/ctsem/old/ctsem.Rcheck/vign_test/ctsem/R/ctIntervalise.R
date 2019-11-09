#' Converts absolute times to intervals for wide format ctsem panel data
#' @param datawide Wide format data, containing absolute time measurements, 
#' to convert to interval time scale. Otherwise as used in \code{\link{ctFit}}.  
#' See \code{\link{ctLongToWide}} to easily convert long format data.
#' @param Tpoints Maximum number of discrete time points (waves of data, or measurement occasions) 
#' for an individual in the input data structure.
#' @param n.manifest number of manifest variables per time point in the data.
#' @param n.TDpred number of time dependent predictors in the data structure.
#' @param n.TIpred number of time independent predictors in the data structure.
#' @param imputedefs if TRUE, impute time intervals based on the measurement occasion (i.e. column)
#' they are in, if FALSE (default), set related observations to NA.  
#' FALSE is recommended unless you are certain that the imputed value 
#' (mean of the relevant time column) is appropriate.  
#' Noise and bias in estimates will result if wrongly set to TRUE.
#' @param manifestNames vector of character strings giving variable names of manifest 
#' indicator variables (without _Tx suffix for measurement occasion).
#' @param TDpredNames vector of character strings giving variable names of time 
#' dependent predictor variables (without _Tx suffix for measurement occasion).
#' @param TIpredNames vector of character strings giving variable names of time 
#' independent predictor variables.
#' @param digits How many digits to round to for interval calculations.
#' @param mininterval set to lower than any possible observed measurement interval, 
#' but above 0 - this is used for filling NA values where necessary and has no 
#' impact on estimates when set in the correct range.  
#' (If all observed intervals are greater than 1, mininterval=1 may be a good choice)
#' @param individualRelativeTime if TRUE (default), the first measurement for each individual is 
#' assumed to be taken at time 0, and all other times are adjusted accordingly.  
#' If FALSE, new columns for an initial wave are created, consisting only of observations 
#' which occurred at the earliest observation time of the entire sample.
#' @param startoffset if 0 (default) uses earliest observation as start time.  
#' If greater than 0, all first observations are NA, with distance of 
#' startoffset to first recorded observation.
#' @details Time column must be numeric!
#' @examples
#' #First load the long format data with absolute times
#' data('longexample')
#' 
#' #Then convert to wide format
#' wideexample <- ctLongToWide(datalong = longexample, id = "id", 
#' time = "time", manifestNames = c("Y1", "Y2", "Y3"), 
#' TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
#' 
#' #Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
#' wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3, 
#' n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
#' TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )
#'  
#' print(wide)
#' @export

ctIntervalise<-function(datawide,Tpoints,n.manifest,n.TDpred=0,n.TIpred=0,imputedefs=F,
  manifestNames='auto', TDpredNames='auto',TIpredNames='auto',
  digits=5,mininterval=.001,individualRelativeTime=TRUE,startoffset=0){

 
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  

  if(n.TDpred > 0 | any(TDpredNames != 'auto')){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  
  if(n.TIpred > 0 | any(TIpredNames != 'auto')){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  tempwide<-as.matrix(datawide,nrow=nrow(datawide)) #ensure matrix for transformations
  
  timeindex<-((Tpoints*n.manifest+(Tpoints)*n.TDpred)+1) : 
    ((Tpoints*n.manifest+(Tpoints)*n.TDpred)+Tpoints) #set appropriate time index for tempwide object
  
  
  if(imputedefs==TRUE) { #impute def vars from column mean
    
    if(any(is.na(colSums(tempwide[,timeindex])))) stop("Time column with no data exists, cannot impute intervals")
    
    message("Warning: imputedefs=T does not make sense (may increase noise and bias estimates) unless the mean of the column times is particularly plausible - use with caution")
    
    tempwide[,timeindex]<-  apply(tempwide[,timeindex],2,function(x) { #set any missing values on time columns to mean of time column
      x[which(is.na(x))]<- mean(x,na.rm=T)
      return(x)})
  }
  
  
  if(imputedefs==FALSE) { #set missing variable values to NA
    message("imputedefs==FALSE (default, recommended) so setting observations with no time value to NA")
    
    if (Tpoints > 1) for(i in 1:(Tpoints-1)){ #for every time, 
      tempwide[which(is.na(tempwide[, timeindex[i]])), #rows that contain missings on time i
        seq(i, #columns that begin with first tdpred at time i 
          n.manifest*Tpoints, #up to last manifest
          Tpoints)] <-  #by tpoints
        NA #if NA, set corresponding manifests to NA
      
      if(n.TDpred>0) {
        tempwide[which(is.na(tempwide[, timeindex[i]])), #rows that contain missings on time i
          seq(n.manifest*Tpoints+i, #columns that begin with first tdpred at time i 
            n.manifest*Tpoints+n.TDpred*Tpoints, #up to last tdpred
            Tpoints)] <-  #by tpoints
          NA #if NA, set corresponding TDpreds to NA
      }
    }
    
    tempwide[which(is.na( #in rows of tempwide where
      tempwide[,timeindex[1]] #the first time points are NA
    )),
      timeindex[1]] <-  mininterval  #set first time point data to mininterval ( we can do this now because we've set variable values to NA already)
    
    
    if (Tpoints > 1) for(i in 2:(Tpoints)){ #for every time after first
      if(any(is.na(tempwide[,timeindex[i]]))){ #if any timing data at time i is missing,
        tempwide[
          which(is.na(tempwide[,timeindex[i]])),
          timeindex[i]
          ] <- #if NA, set time to mininterval + earlier time (has no observation so time can be arbitrarily set as long as it is between neighbouring observations)
          mininterval + tempwide[
            which(is.na(tempwide[,timeindex[i]])),
            timeindex[i]-1
            ] 
      }
    }
  }  
  
  
  
  
  
  if(nrow(tempwide)==1) individualRelativeTime<-TRUE #if only 1 subject then flatten start time to 0
  if(length(unique(tempwide[,Tpoints*n.manifest+(Tpoints)*n.TDpred+1])) == 1) individualRelativeTime<-TRUE #if all Tpoint 1 times are equal, then flatten start time to 0
  
  if(n.TIpred>0) tempwide<-tempwide[,-ncol(tempwide) : 
      -(ncol(tempwide)-n.TIpred+1),drop=FALSE] #remove TI predictors for moment
  
  if(all(is.na(
    tempwide[,Tpoints*n.manifest+(Tpoints)*n.TDpred+1])
  )){ #if no data in first time column
    
    message("first time column empty! setting to startoffset and adjusting")
    
    tempwide[,Tpoints*n.manifest+(Tpoints)*n.TDpred+1] <- startoffset
    individualRelativeTime<-TRUE
  }
  
  
  
  
  
  if(individualRelativeTime==FALSE){ #if there are multiple cases and the observations do not all start at the same time, and the wave structure should be retained
    
    
    temp <- matrix(c( rep(NA, times = (n.manifest) * nrow(tempwide)) , #add blank first observation columns to front 
      tempwide[,1:(n.manifest*Tpoints)], #then manifests
      rep(NA, times = (n.TDpred) * nrow(tempwide)), #then extra blank TDpredictors columns as needed
      tempwide[,(n.manifest*Tpoints+1) : 
          (n.manifest*Tpoints+n.TDpred*(Tpoints)+Tpoints)] #then rest of data
    ), nrow=nrow(tempwide))
    
    Tpoints<-Tpoints+1 #because extra column was added
    message('Extra measurement occasion created in data structure -- Tpoints now ', Tpoints)
    
    colnames(temp)<- ctWideNames(n.manifest=n.manifest, n.TDpred = n.TDpred,
      Tpoints=Tpoints, manifestNames=manifestNames, TDpredNames=TDpredNames, TIpredNames=TIpredNames, n.TIpred=0) #set colnames here for easier debugging, but set later too
    
    timeindex<-((Tpoints*n.manifest+(Tpoints)*n.TDpred)+1) : 
      ((Tpoints*n.manifest+(Tpoints)*n.TDpred)+Tpoints-1) #set index of time variables in temp object
    
    if(startoffset>0) temp[,timeindex] <- temp[,timeindex]+startoffset   #if setting start offset, set first interval to start offset (first obs will then all be NA)
    
    
    if(startoffset==0){ #if we want the first column block (rather than the second) to contain the earliest observations
      
      temp[,timeindex] <- temp[,timeindex]-
        temp[which(temp[,timeindex[1]] == 
            min(temp[,timeindex[1]],na.rm=T)),timeindex[1]][1] #subtract min time from all intervals to 0 beginning
      
      temp[which(temp[,timeindex[1]]==0),1:n.manifest] <-  temp[which(temp[,timeindex[1]]==0),
        (n.manifest+1):(n.manifest*2)] #any manifest at time 0 go to first column block
      temp[which(temp[,timeindex[1]]==0),(n.manifest+1):(n.manifest*2)]<-NA #and corresponding manifest in second column block set to NA
      
      if(n.TDpred>0) {
        temp[which(temp[,timeindex[1]]==0),(n.manifest*Tpoints+1) : 
            (n.manifest*Tpoints+n.TDpred)] <- #any TDpreds at time 0 go to first TDpred column
          temp[which(temp[,timeindex[1]]==0),
            (n.manifest*Tpoints+1+n.TDpred) : 
              (n.manifest*Tpoints+n.TDpred*2)]
        
        temp[which(temp[,timeindex[1]]==0),
          (n.manifest*Tpoints+1+n.TDpred) : 
            (n.manifest*Tpoints+n.TDpred*2)] <- NA #with corresponding TDpreds in 2nd block set NA
      }
      
      temp[which(temp[,timeindex[1]]==0),
        timeindex[1]] <- mininterval #and corresponding interval for second column block set to mininterval
    }
  }
  
  
  
  
  
  
  if(individualRelativeTime==TRUE){ #if the user wants to set the first obs time for each subject to 0 (generally makes sense - only problematic if some effect directly on the wave is to be implemented)
    
    timeindex <- ((Tpoints*n.manifest+(Tpoints)*n.TDpred)+1) : 
      (Tpoints*n.manifest+(Tpoints)*n.TDpred+Tpoints) #set appropriate time index for temp object
    temp <- tempwide
    
    if(any(is.na(temp[,timeindex[1]]))) temp[is.na(temp[,timeindex[1]]),timeindex[1]]<-0 #Any missing time data for first obs are set to 0 - this is ok to set because we've removed data or imputed missing time obs above
    intervals<-temp[,timeindex,drop=FALSE] #extract intervals for a moment so structure is not broken
    
    intervals<-matrix(apply(intervals,1,function(x) {
      x<-c(x)#subtract first obs time from all obs to flatten start time
      x <- x - as.numeric(x[1])
      return(x)
    })
      ,ncol=(Tpoints),byrow=TRUE)
    
    temp[,timeindex]<-intervals #push intervals back into data structure
    temp<-temp[,-timeindex[1],drop=FALSE] #remove first time column (as this is now 0 in all cases)
    timeindex<-timeindex[-length(timeindex)] #adjust timeindex accordingly
  }
  
  
  #adjust absolute times to represent intervals
  if (Tpoints > 2) for(i in (Tpoints-1):2) { #for every time obs from the last to the 2nd
    temp[,timeindex[i]] <-  temp[,timeindex[i]] - temp[,timeindex[i-1]] #set the obs to the difference between itself and the next earlier time
  }
  
  temp[,timeindex] <- round(temp[,timeindex],digits=digits) #round any intervals to specified digits
  
  #   temp[apply(temp[,timeindex],1,function(x) all(is.na(x))),timeindex]<- mininterval #set any rows with all missing times to intervals of mininterval
  
  
  if(n.TIpred>0) temp <- cbind(temp, datawide[,(ncol(datawide)-n.TIpred+1) :  ncol(datawide),drop=FALSE]) #add TI predictors back  
  

  colnames(temp)<-ctWideNames(n.manifest=n.manifest, n.TDpred = n.TDpred,
    Tpoints=Tpoints, manifestNames=manifestNames,TDpredNames=TDpredNames,TIpredNames=TIpredNames, n.TIpred=n.TIpred)
  
  return(temp)
}
