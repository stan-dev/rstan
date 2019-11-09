#' Fits a multiple group continuous time model.
#' 
#' Fits a single continuous time structural equation models to multiple groups (where each group contains 1 or more subjects),
#' by default, all parameters are free across groups.  Can also be used to easily estimate seperate models for each group.
#' 
#' @param dat Wide format data, as used in \code{\link{ctFit}}.  See \code{\link{ctLongToWide}} to
#' easily convert long format data.
#' @param groupings For wide format: Vector of character labels designating group membership for each row of dat.  
#' For long format: Named list of groups, with each list element containing a vector of subject id's for the group.
#' In both cases, group names will be prefixed on relevant parameter estimates in the summary.
#' @param ctmodelobj Continuous time model to fit, specified via \code{\link{ctModel}} function.
#' @param dataform either "wide" or "long" depending on which input format you wish to use for the data. 
#' See details of \code{\link{ctFit}} and or vignette.
#' @param fixedmodel Modified version of ctmodelobj, wherein any parameters you wish to keep 
#' fixed over groups should be given the value 'groupfixed'.  
#' If specified, all other parameters will be free across groups.
#' @param freemodel Modified version of ctmodelobj, wherein any parameters you wish to free across groups
#' should be given the label 'groupfree'.  
#' If specified, all other parameters will be fixed across groups.  
#' If left NULL, the default, all parameters are free across groups.
#' @param showInits Displays start values prior to optimization
#' @param omxStartValues A named vector containing the raw (potentially log transformed) OpenMx starting values for free parameters, as captured by
#' OpenMx function \code{omxGetParameters(ctmodelobj$mxobj)}. These values will take precedence 
#' over any starting values already specified using ctModel.
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to discourage parameters from boundary conditions, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help / speed optimization, though results in user specified inits being ignored for the final fit.
#' @param retryattempts Number of fit retries to make.
#' @param ... additional arguments to pass to \code{\link{ctFit}}.
#' @return Returns an OpenMx fit object.
#' @details Additional \code{\link{ctFit}} parameters may be specified as required. Confidence intervals for any matrices and or parameters 
#' may be estimated afer fitting using \code{\link{ctCI}}.
#' 
#' @examples 
#' \dontrun{
#' 
#' #Two group model, all parameters except LAMBDA[3,1] constrained across groups.
#' data(ctExample4)
#' basemodel<-ctModel(n.latent=1, n.manifest=3, Tpoints=20,
#'                    LAMBDA=matrix(c(1, 'lambda2', 'lambda3'), nrow=3, ncol=1),
#'                    MANIFESTMEANS=matrix(c(0, 'manifestmean2', 'manifestmean3'), 
#'                    nrow=3, ncol=1), TRAITVAR = 'auto')
#' 
#' freemodel<-basemodel
#' freemodel$LAMBDA[3,1]<-'groupfree'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(dat=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, freemodel=freemodel)
#' summary(multif,group=1)
#' 
#' 
#' 
#' #fixed model approach
#' fixedmodel<-basemodel
#' fixedmodel$LAMBDA[2,1]<-'groupfixed'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(dat=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, fixedmodel=fixedmodel)
#' summary(multif) 
#'}
#' 
#' 
#' @seealso \code{\link{ctFit}} and \code{\link{ctModel}}
#' @export


ctMultigroupFit<-function(dat,groupings,ctmodelobj,dataform='wide',fixedmodel=NA,freemodel=NA,
 carefulFit=TRUE, omxStartValues=NULL,
  retryattempts=15,showInits=FALSE,...){

  if(dataform=='wide') if(any(suppressWarnings(!is.na(as.numeric(groupings))))) stop("grouping variable must not contain purely numeric items")
  if(dataform=='wide') if(length(groupings)!= nrow(dat)) stop('length of groupings does not equal number of rows of dat')
  if(dataform=='long') if(length(unlist(groupings)) != length(unique(dat[,'id']))) stop('groupings list does not contain the right number of subjects!')
  if(dataform=='long') if(any(suppressWarnings(!is.na(as.numeric(names(groupings)))))) stop("grouping variable must not contain purely numeric items")
  
  if(all(is.na(fixedmodel))) fixedmodel<-ctmodelobj #so that it is not null or na
  
  startparams<-c() #to fill as needed
  
  
  
  omxmodels<-list() #blank list preparing for model input
  if(dataform=='wide') groups <- unique(groupings)
  if(dataform=='long') groups <- unique(names(groupings))
  for(i in groups){ #for every specified group

    if(dataform=='wide') singlegroup<-dat[which(groupings == i),,drop=FALSE] #data for the group
    if(dataform=='long') singlegroup<-dat[which(dat[,'id'] %in% groupings[[i]]),,drop=FALSE] #data for the group
    
    singlectspec<-ctmodelobj
    
    if(all(is.na(freemodel))) freemodel <- lapply(ctmodelobj,function(x) { x<-rep('groupfree',length(x))}) #if no freemodel specified then free all params at this point
    
    for(m in 1:length(ctmodelobj)) { #for every element of the ctmodelobj list
      if(is.matrix(ctmodelobj[[m]])){ #if the element is a matrix
      for(j in 1:length(ctmodelobj[[m]])){ #for every slot in the matrix
        if(freemodel[[m]][j]=="groupfree"){ #if the slot is free in freemodel and not fixed in fixedmodel
        jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
        if(!is.na(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither NA or fixed to a value, then
          singlectspec[[m]][j] <- paste0(i,'_',ctmodelobj[[m]][j]) #give the label a group specific prefix
          
          
        }
      }
        if(any(!is.na(fixedmodel))){
        if(fixedmodel[[m]][j]=="groupfixed"){
          jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
          if(!is.null(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither null or fixed to a value, then
            singlectspec[[m]][j] <- paste0(ctmodelobj[[m]][j]) #use the global label
          }
        }
      }
    }
      }
    }
   

    
    
    
    if(carefulFit==TRUE) message('Begin carefulFit start value estimation for group ', i)
    
    omxmodel<-ctFit(singlegroup,singlectspec,nofit=TRUE, carefulFit=carefulFit,dataform=dataform,...) #omxmodel for group i
    ctfitargs<-omxmodel$ctfitargs
    omxmodel<-omxmodel$mxobj
    
    if(carefulFit==TRUE) {
      startparams<-c( startparams[ !( names(startparams) %in%  #get inits
          names(OpenMx::omxGetParameters(omxmodel))) ], #that are not found in the new fits inits
        OpenMx::omxGetParameters(omxmodel) ) #and combine the two vectors
    }
    
    omxmodel<- OpenMx::mxRename(omxmodel, newname=i) #change name of omxmodel for group i
      
  omxmodels[[i]]<-omxmodel #if fitting single multigroup model, add omxmodel for group i to list of omxmodels for all groups
    
  } #end loop over groups
  
  
  
    
    fullmodel <- OpenMx::mxModel('ctsem multigroup', #output multigroup omxmodel
      mxFitFunctionMultigroup(c(paste0(groups))),
#       mxComputeSequence(list(
#         mxComputeGradientDescent(gradientAlgo="central", nudgeZeroStarts=FALSE, 
#           maxMajorIter=1000, gradientIterations = 1),
        # mxComputeReportDeriv(),
      omxmodels)
    

    fullmodel<-OpenMx::omxAssignFirstParameters(fullmodel)

    if(!is.null(omxStartValues)) fullmodel<-omxSetParameters(fullmodel,
      labels=names(omxStartValues)[names(omxStartValues) %in% names(omxGetParameters(fullmodel))],
      values=omxStartValues[names(omxStartValues) %in% names(omxGetParameters(fullmodel))],strict=FALSE)
    
      fullmodel<-OpenMx::mxTryHard(fullmodel,initialTolerance=1e-14,
      showInits=showInits,
      bestInitsOutput=FALSE,
      extraTries=retryattempts,loc=1,scale=.2,paste=FALSE) 
    
      fullmodel<-list(mxobj=fullmodel, ctfitargs=ctfitargs, ctmodelobj=ctmodelobj, groups=groups)
      class(fullmodel)<-'ctsemMultigroupFit'
      
      
    return(fullmodel)
  
}


