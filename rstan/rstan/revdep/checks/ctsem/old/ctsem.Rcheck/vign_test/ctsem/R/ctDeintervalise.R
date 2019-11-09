#' ctDeintervalise
#' 
#' Converts intervals in ctsem long format data to absolute time
#' @param datalong data to use, in ctsem long format (attained via function ctWideToLong)
#' @param id character string denoting column of data containing numeric identifier for each subject.
#' @param dT character string denoting column of data containing time interval preceding observations in that row.
#' @param startoffset Number of units of time to offset by when converting.
#' @export
ctDeintervalise<-function(datalong,id='id', dT='dT',startoffset=0){

  message(paste0("Converting intervals to absolute time:  Any missing intervals on 1st row of each subject are assumed to occur at earliest measurement time (", startoffset ,"), any other missing intervals render subsequent intervals for the subject unusable so time variables are set NA"))

  # initialmissingcount <- ifelse(is.na(datalong[1,dT]),1,0)
  # othermissingcount<-0
  datalong[1,dT]<-sum(c(datalong[1,dT],startoffset),na.rm=TRUE) #datalong row 1 equals first interval and offset
  
  for(i in 2:nrow(datalong)){ #for subsequent rows
          if(datalong[i,id]==datalong[i-1,id]){ #check if the subject is the same as the row above
            # othermissingcount <- ifelse(is.na(datalong[i,id]==datalong[i-1,id]),othermissingcount+1, othermissingcount)            
            datalong[i,dT]<-sum(datalong[(i-1):i,dT],na.rm=FALSE) #if same subject, sum the new interval with the prev total time
      } else {
        # initialmissingcount <- ifelse(is.na(datalong[i,dT]),initialmissingcount+1,initialmissingcount)
        datalong[i,dT]<-sum(c(datalong[i,dT],startoffset),na.rm=T) #otherwise create new total time with new interval and offset
      }
  }
  colnames(datalong)[colnames(datalong) %in% dT] <-'time'
  return(datalong)
}
