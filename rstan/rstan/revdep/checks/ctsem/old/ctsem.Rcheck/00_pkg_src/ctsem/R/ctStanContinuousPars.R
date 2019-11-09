#'ctStanContinuousPars
#'
#'Returns the continuous time parameter matrices for specified subjects of a ctStanFit fit object
#'
#'@param ctstanfitobj fit object from \code{\link{ctStanFit}}
#'@param subjects Either 'all', or integers denoting which subjects to perform the calculation over. 
#'When multiple subjects are specified, the returned matrices will be a mean over subjects.
#'@param iter Either character string 'all' which will then use all post-warmup iterations, 
#'or an integer specifying which iteration/s to use.
#'@param calcfunc Function to apply over samples, must return a single value. 
#'By default the median over all samples is returned using the \code{\link[stats]{quantile}} function, 
#'but one might also be interested in the \code{\link[base]{mean}} or \code{\link[stats]{sd}}, for instance.
#'@param calcfuncargs A list of additional parameters to pass to calcfunc. 
#'For instance, with the default of calcfunc = quantile, 
#'the probs argument is needed to ensure only a single value is returned.
#'@examples
#'\dontrun{
#'#posterior median over all subjects (also reflects mean of unconstrained pars)
#'ctStanContinuousPars(ctstantestfit)
#'
#'#posterior 97.5% quantiles for subject 2
#'ctStanContinuousPars(ctstantestfit, subjects=2, calcfunc=quantile, 
#'calcfuncargs=list(probs=0.975))
#'}
#'@export
ctStanContinuousPars <- function(ctstanfitobj,subjects='all',iter='all',
  calcfunc=quantile,calcfuncargs=list(probs=0.5)){
  if(subjects[1] != 'all' && !is.integer(as.integer(subjects))) stop('
    subjects argument must be either "all" or an integer denoting specific subjects')
  
  if(class(ctstanfitobj)!='ctStanFit') stop('Not an object of class ctStanFit')
  
  e<-extract(ctstanfitobj) #first dim of subobjects is iter, 2nd subjects
  niter=dim(e$DRIFT)[1]
  
  if(iter!='all') {
    if(!any(iter %in% 1:niter)) stop('Invalid iteration specified!')
    e=lapply(e, function(x) {
      xdims=dim(x)
      out=array(eval(parse(text=
          paste0('x[iter',if(length(xdims)>1) paste0(rep(',',length(xdims)-1),collapse=''),']')
      )),dim=c(length(iter),xdims[-1]))
      return(out)
    }
    )
  }
  
  nsubjects <- ctstanfitobj$data$nsubjects #dim(e$indparams)[2]
  # if(is.null(nsubjects)) nsubjects=1
  
  # if(subjects[1]=='all') subjects=1:nsubjects
  
  if(subjects=='all') collapsemargin <- 1 else collapsemargin<-c(1,2)
  # if(collapseIterations) collapsemargin=1
  
  for(matname in c('DRIFT','DIFFUSION','asymDIFFUSION','CINT','T0MEANS', 
    'T0VAR','MANIFESTMEANS',if(!is.null(e$pop_MANIFESTVAR)) 'MANIFESTVAR','LAMBDA', if(!is.null(e$pop_TDPREDEFFECT)) 'TDPREDEFFECT')){

    subindex <- ctstanfitobj$data[[paste0(matname,'subindex')]]
    # if(dim(e[[matname]])[2] > 1) subselection <- subjects else subselection <- 1
    if(max(subindex) > 1 || subjects =='all') subselection <- 1:max(subindex) else subselection <- 1
    subselection=ctstanfitobj$standata[[paste0(matname,'subindex')]][subselection]
    
    # vector <- FALSE
    
    # if(matname %in% c('T0MEANS','CINT', 'MANIFESTMEANS')) vector <- TRUE
    
    calcfuncargs$collapsemargin = collapsemargin
    calcfuncargs$collapsefunc=calcfunc

    # if(!vector) {
      calcfuncargs$inarray = e[[ifelse(subjects!='all', matname, paste0('pop_',matname))]]
      if(subjects!='all') calcfuncargs$inarray  <- calcfuncargs$inarray[,subselection,,,drop=FALSE]
      assign(matname, array(do.call(ctCollapse,calcfuncargs),
        dim=dim(e[[ifelse(subjects!='all', matname, paste0('pop_',matname))]])[-c(1,if(subjects!='all') 2)]))
    # }
    
    # if(vector) {
    #   calcfuncargs$inarray = e[[matname]][,subselection,,drop=FALSE]
    #   assign(matname, array(do.call(ctCollapse,calcfuncargs),dim=c(dim(e[[matname]])[-c(1,2)],1)))
    # }
    
  }
  
  # DRIFTHATCH<-DRIFT %x% diag(nrow(DRIFT)) + diag(nrow(DRIFT)) %x% DRIFT
  # asymDIFFUSION<-matrix(-solve(DRIFTHATCH) %*% c(DIFFUSION), nrow=nrow(DRIFT)) 
  
  ln=ctstanfitobj$ctstanmodel$latentNames
  mn=ctstanfitobj$ctstanmodel$manifestNames
  tdn=ctstanfitobj$ctstanmodel$TDpredNames
  dimnames(DRIFT)=list(ln,ln)
  dimnames(DIFFUSION)=list(ln,ln)
  dimnames(asymDIFFUSION)=list(ln,ln)
  rownames(CINT)=ln
  rownames(MANIFESTMEANS)=mn
  rownames(T0MEANS)=ln
  
  dimnames(T0VAR)=list(ln,ln)
  dimnames(asymDIFFUSION)=list(ln,ln)
  dimnames(LAMBDA)=list(mn,ln)
  
  model<-list(DRIFT=DRIFT,T0VAR=T0VAR,DIFFUSION=DIFFUSION,asymDIFFUSION=asymDIFFUSION,CINT=CINT,T0MEANS=T0MEANS,
    MANIFESTMEANS=MANIFESTMEANS, LAMBDA=LAMBDA)
  
  if(!is.null(e$pop_MANIFESTVAR)) {
    dimnames(MANIFESTVAR)=list(mn,mn)
    model$MANIFESTVAR=MANIFESTVAR
    
  }
  
  if(!is.null(e$pop_TDPREDEFFECT)) {
    dimnames(TDPREDEFFECT)=list(ln,tdn)
    model$TDPREDEFFECT<-TDPREDEFFECT
  }
  
  return(model)
}

