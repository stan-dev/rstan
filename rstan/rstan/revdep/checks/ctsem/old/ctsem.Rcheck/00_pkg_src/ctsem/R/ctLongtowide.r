#' ctLongToWide
#' Restructures time series / panel data from long format to wide format for ctsem analysis 
#' @param datalong dataset in long format, including subject/id column, observation time 
#' (or change in observation time, with 0 for first observation) column, 
#' indicator (manifest / observed) variables, 
#' any time dependent predictors, and any time independent predictors.
#' @param id character string giving column name of the subject/id column
#' @param time character string giving column name of the time columnn
#' @param manifestNames vector of character strings giving column names of manifest indicator variables
#' @param TDpredNames vector of character strings giving column names of time dependent predictor variables
#' @param TIpredNames vector of character strings giving column names of time independent predictor variables
#' @details Time column must be numeric
#' @export 
#' @seealso \code{\link{ctIntervalise}}
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


ctLongToWide <- function(datalong, id, time, manifestNames, TDpredNames=NULL, TIpredNames=NULL){

  data_long 			<- as.data.frame(datalong)[c(id, 
    manifestNames, 
    TDpredNames, 
    time, 
    TIpredNames)]
  names(data_long) 	<- c("id", 
    manifestNames, 
    TDpredNames, 
    "time",
    TIpredNames)
  
  if(any(is.na(data_long[,'time']))){
    message(paste0('Observations with missing time information found - removing ', sum(is.na(data_long[,'time'])), ' rows'))
    data_long<-data_long[!is.na(data_long[,'time']),]
  }

  data_long <- data_long[order(data_long[,'id'],data_long[,'time']),]  # order by id

  discrete.time.point<-rep(1,nrow(data_long))

  for(i in 2:nrow(data_long)){ #number discrete time points
    if(data_long[i,"id"]==data_long[(i-1),"id"]) {
      discrete.time.point[i] <- discrete.time.point[i-1] + 1
    }
  }
  
  data_long<-cbind(data_long,discrete.time.point) #add discrete time point column  
  
  
  manifestNames_wide<-stats::reshape(data_long,  #create wide format manifestNames
    v.names = manifestNames, 
    idvar = "id", 
    timevar = "discrete.time.point", 
    direction = "wide",
    drop=c(TDpredNames,"time",TIpredNames)) [,-1,drop=FALSE]
  
  if(is.list(manifestNames_wide)) {
    manifestNames_wide<-as.matrix(manifestNames_wide,nrow=1) #because reshape outputs lists when only one row!
    colnames(manifestNames_wide)<-paste0(manifestNames,rep(paste0('_T',0:( (ncol(manifestNames_wide)/length(manifestNames))-1)),each=length(manifestNames)))
  }
  
  
  Tpoints<-ncol(manifestNames_wide)/length(manifestNames)
  message(paste0("Discrete time points created:  ",Tpoints)) #report discrete time points
  
  
  
  TDpred_wide<-matrix(NA,ncol=0,nrow=nrow(manifestNames_wide)) #set null dataframe to avoid errors
  if(!is.null(TDpredNames[1])){ #if TDpredNames are specified
    TDpred_wide<-stats::reshape(data_long,  #create wide format manifestNames
      v.names = TDpredNames, 
      idvar = "id", 
      timevar = "discrete.time.point", 
      direction = "wide",
      drop=c(manifestNames,"time",TIpredNames)) [,-1,drop=FALSE]
    
    if(is.list(TDpred_wide)) {
      TDpred_wide<-as.matrix(TDpred_wide,nrow=1) #because reshape outputs lists when only one row!
      colnames(TDpred_wide)<-paste0(TDpredNames,
        rep(paste0('_T',0:( (ncol(TDpred_wide)/length(TDpredNames))-1)),each=length(TDpredNames)))
    }
    
    
    # for(i in 1:length(TDpredNames)){ #for every TDpredictor
    #   predi <- stats::reshape(data_long, 
    #     v.names = TDpredNames[i], 
    #     idvar = "id", 
    #     timevar = "discrete.time.point", 
    #     direction = "wide",
    #     drop=c(manifestNames,TDpredNames[-i],"time",TIpredNames)) [,-1]
    #   
    #   if(!is.data.frame(predi)) predi<-as.matrix(predi,nrow=1)
    #   colnames(predi) <- paste0(TDpredNames,'_T',0:(ncol(predi)-1))
    #   
    # 
    #   
    #   TDpred_wide<-cbind(TDpred_wide,
    #     predi) #bind TD predictors variable i to full TDpredictor object
    # }  
  }
  
  
  
  time_wide<-stats::reshape(data_long,  
    idvar = "id",  #create time variables block
    timevar = "discrete.time.point", 
    direction = "wide",
    drop=c(manifestNames,TDpredNames,TIpredNames)) [,-1,drop=F]
  if(is.list(time_wide)) time_wide<-as.matrix(time_wide,nrow=1) #because reshape outputs lists when only one row!
  colnames(time_wide) <- paste0('T',0:(ncol(time_wide)-1))

  
  
  
  
  TIpred_wide<-data.frame(matrix(NA,ncol=length(TIpredNames),nrow=nrow(manifestNames_wide))) 
  if(!is.null(TIpredNames[1])){ #if there are time independent predictors
    message('Extracting first non-missing value for time independent predictors')
    # pb<-utils::txtProgressBar(1,nrow(data_long))

    idi<- -1
    widei<-0
    for(longi in 1:nrow(data_long))  {
      if(data_long[longi,'id'] != idi){ #if new id
        idi<-data_long[longi,'id'] #update id
        widei<-widei+1 #increment wide row count
        TIpred_wide[widei,] <- data_long[longi,TIpredNames] #update wide with long info
      }
      if(data_long[longi,'id'] == idi && any(is.na(TIpred_wide[widei,]))){ #if same id and some wide values NA, try this long row
        TIpred_wide[widei,is.na(TIpred_wide[widei,])]<-data_long[longi,TIpredNames[is.na(TIpred_wide[widei,])]]
      }
        # utils::setTxtProgressBar(pb,longi)  
    }
       
      
      
      
 colnames(TIpred_wide)<-TIpredNames
  }
  
  
  
  data_wide1<-as.matrix(cbind(manifestNames_wide,
    TDpred_wide,
    time_wide,
    TIpred_wide)) #combine data types together for ctsem wide format 
  
  rownamevector<- unique(data_long$id) #set row names to id's
  dimnames(data_wide1)[[1]]<-rownamevector #and put the row names onto the object (errors generated?)
  
  return(data_wide1)
}
