#' Fit a ctsem object
#' 
#' This function fits continuous time SEM models specified via \code{\link{ctModel}}
#' to a dataset containing one or more subjects.
#' @param dat the data you wish to fit a ctsem model to, in either wide format (one individual per row), 
#' or long format (one time point of one individual per row). See details. 
#' @param dataform either "wide" or "long" depending on which input format you wish to use for the data. See details and or vignette.
#' @param ctmodelobj the ctsem model object you wish to use, specified via the \code{\link{ctModel}} function.
#' @param fit if FALSE, output only openmx model without fitting
#' @param nofit Deprecated. If TRUE, output only openmx model without fitting
#' @param objective 'auto' selects either 'Kalman', if fitting to single subject data, 
#' or 'mxRAM' for multiple subjects. For single subject data, 'Kalman' uses the \code{mxExpectationStateSpace }
#' function from OpenMx to implement the Kalman filter. 
#' For more than one subject, 'mxRAM' specifies a wide format SEM with a row of data per subject.
#' 'cov' may be specified, in which case the 'meanIntervals' argument is set to TRUE, and the covariance matrix
#' of the supplied data is calculated and fit instead of the raw data. This is much faster but only a rough approximation,
#' unless there are no individual differences in time interval and no missing data.
#' 'Kalman' may be specified for multiple subjects, however as no trait matrices are used by the Kalman filter
#' one must consider how average level differences between subjects are accounted for.
#' See \code{\link{ctMultigroupFit}} for the possibility to apply the Kalman filter over multiple subjects)
#' @param stationary Character vector of T0 matrix names in which to constrain any 
#' free parameters to stationarity. 
#' Defaults to \code{c('T0TRAITEFFECT','T0TIPREDEFFECT')}, constraining only
#' between person effects to stationarity. Use \code{NULL} for no constraints,
#' or 'all' to constrain all T0 matrices.
#' @param optimizer character string, defaults to the open-source 'CSOLNP' optimizer that is distributed
#' in all versions of OpenMx. However, 'NPSOL' may sometimes perform better for these problems,
#' though requires that you have installed OpenMx via the OpenMx web site, by running:
#' \code{source('http://openmx.psyc.virginia.edu/getOpenMx.R')} 
#' @param showInits if TRUE, prints the list of 
#' starting values for free parameters. These are the 'raw' values used by OpenMx, 
#' and reflect the log (var / cov matrices) or -log(DRIFT matrices) transformations used in ctsem.
#' These are saved in the fit object under \code{fitobject$omxStartValues}.
#' @param retryattempts Number of times to retry the start value randomisation and fit procedure, if non-convergance or uncertain fits occur.
#' @param iterationSummary if TRUE, outputs limited fit details after every fit attempt.
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to force MANIFESTVAR, DRIFT, TRAITVAR, MANIFESTTRAITVAR parameters to remain close to 0, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help to ensure optimization begins at sensible, non-exteme values, 
#' though results in any user specified start values being ignored for the final fit (though they are still used for initial fit).
#' @param carefulFitWeight Positive numeric. Sets the weight for the penalisation (or prior) applied by the carefulFit algorithm. 
#' Generally unnecessary to adjust, may be helpful to try a selection of values (perhaps between 0 and 1000) when optimization is problematic.
#' @param plotOptimization If TRUE, uses checkpointing for OpenMx function \code{mxRun}, set to checkpoint every iteration, 
#' output checkpoint file to working directory, then creates a plot for each parameter's values over iterations.
#' @param meanIntervals Use average time intervals for each column for calculation 
#' (both faster and inaccurate to the extent that intervals vary across individuals).
#' @param discreteTime Estimate a discrete time model - ignores timing information, parameter
#' estimates will correspond to those of classical vector autoregression models, 
#' OpenMx fit object will be directly output, thus ctsem summary and plot functionality will be unavailable.
#' Time dependent predictor type also becomes irrelevant.
#' @param asymptotes when TRUE, optimizes over asymptotic parameter matrices instead of continuous time parameter matrices. 
#' Can be faster for optimization and in some cases makes reliable convergance easier. Will result in equivalent models 
#' when continuous time input matrices (DRIFT, DIFFUSION, CINT) are free, but fixing the values of 
#' any such matrices will result in large differences - a value of 0 in a cell of the normal continuous time DIFFUSION matrix
#' does not necessarily result in a value of 0 for the asymptotic DIFFUSION matrix, for instance.
#' @param verbose Integer between 0 and 3. Sets mxComputeGradientDescent messaging level, defaults to 0.
#' @param useOptimizer Logical. Defaults to TRUE.  Passes argument to \code{mxRun}, 
#' useful for using custom optimizers or fitting to specified parameters.
#' @param omxStartValues A named vector containing the raw (potentially log transformed) OpenMx starting values for free parameters, as captured by
#' OpenMx function \code{omxGetParameters(ctmodelobj$mxobj)}. These values will take precedence 
#' over any starting values already specified using ctModel.
#' @param transformedParams Logical indicating whether or not to log transform 
#' certain parameters internally to allow unconstrained estimation over
#' entire 'sensible' range for parameters. 
#' When TRUE (default) raw OpenMx parameters (only reported if \code{verbose=TRUE} argument used
#' for summary function) will reflect these transformations and may be harder to 
#' interpret, but summary matrices are reported as normal.
#' @param crossEffectNegStarts Logical. If TRUE (default) free DRIFT matrix cross effect parameters have starting values 
#' set to small negative values (e.g. -.05), if FALSE, the start values are 0. The TRUE setting is useful for easy 
#' initialisation of higher order models, while the FALSE setting is useful when one has already estimated a model without cross effects,
#' and wishes to begin optimization from those values by using the omxStartValues switch.
#' are re-transformed into regular continuous time parameter matrices, and may be interpreted as normal.
#' @param datawide included for compatibility with scripts written for earlier versions of ctsem. 
#' Do not use this argument, instead use the dat argument, and the dataform argument now specifies whether the
#' data is in wide or long format.
#' @details For full discussion of how to structure the data and use this function, see the vignette using: \code{vignette('ctsem')}, or
#' the data examples \code{data("longexample") ; longexample} for long and \code{data("datastructure") ; datastructure} for wide. 
#' If using long format, the subject id column must be numeric and grouped by ascending time within subject, and named 'id'. 
#' The time column must also be numeric, and representing absolute time (e.g., since beginning of study, *not* time intervals),
#' and called 'time'.
#' Models are specified using the \code{\link{ctModel}} function.
#' For help regarding the summary function, see \code{\link{summary.ctsemFit}}, 
#' and for the plot function, \code{\link{plot.ctsemFit}}.
#' Multigroup models may be specified using \code{\link{ctMultigroupFit}}.
#' Confidence intervals for any matrices and or parameters 
#' may be estimated using \code{\link{ctCI}}.
#' Difficulties during estimation can sometimes be alleviated using \code{\link{ctRefineTo}} instead of \code{\link{ctFit}} -- 
#' this uses a multistep fit procedure.
#' @examples 
#' ## Examples set to 'dontrun' because they take longer than 5s.
#' \dontrun{
#' mfrowOld<-par()$mfrow
#' par(mfrow=c(2, 3))
#' 
#' ### example from Driver, Oud, Voelkle (2017), 
#' ### simulated happiness and leisure time with unobserved heterogeneity.
#' data(ctExample1)
#' traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
#'   manifestNames=c('LeisureTime', 'Happiness'), 
#'   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' traitfit <- ctFit(dat=ctExample1, ctmodelobj=traitmodel)
#' summary(traitfit)
#' plot(traitfit, wait=FALSE)
#' 
#' 
#' ###Example from Voelkle, Oud, Davidov, and Schmidt (2012) - anomia and authoritarianism.  
#' data(AnomAuth) 
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#' Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#' summary(AnomAuthfit)
#' 
#' 
#' ### Single subject time series - using Kalman filter (OpenMx statespace expectation)
#' data('ctExample3')
#' model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
#'   LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1), 
#'   CINT= matrix('cint'),
#'   MANIFESTMEANS = matrix(c(0, 'manifestmean2', 'manifestmean3'), nrow = 3, 
#'     ncol = 1))
#' fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = 'Kalman', 
#'   stationary = c('T0VAR'))
#' 
#' 
#' ###Oscillating model from Voelkle & Oud (2013). 
#' data("Oscillating")
#' 
#' inits <- c(-39, -.3, 1.01, 10.01, .1, 10.01, 0.05, .9, 0)
#' names(inits) <- c("crosseffect","autoeffect", "diffusion",
#'   "T0var11", "T0var21", "T0var22","m1", "m2", 'manifestmean')
#' 
#' oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11,
#'   MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1),
#'   LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
#'   T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1),
#'   T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
#'   DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2),
#'   CINT = matrix(0, ncol = 1, nrow = 2),
#'   MANIFESTMEANS = matrix('manifestmean', nrow = 1, ncol = 1),
#'   DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2),
#'   startValues=inits)
#' 
#' oscillatingf <- ctFit(Oscillating, oscillatingm, carefulFit = FALSE)
#' }
#' @import OpenMx
#' @export

ctFit  <- function(dat, ctmodelobj, dataform='wide',
  objective='auto', 
  stationary=c('T0TRAITEFFECT','T0TIPREDEFFECT'), 
  optimizer='CSOLNP', 
  retryattempts=15, iterationSummary=FALSE, carefulFit=TRUE,  
  carefulFitWeight=100,
  showInits=FALSE, asymptotes=FALSE,
  meanIntervals=FALSE, plotOptimization=FALSE, 
  crossEffectNegStarts=TRUE,
  fit = TRUE, nofit=FALSE, discreteTime=FALSE, verbose=0, useOptimizer=TRUE,
  omxStartValues=NULL, transformedParams=TRUE,
  datawide=NA){
  
  if(length(datawide) > 1 || !is.na(datawide)) {
    warning('ctsem now uses the dat argument instead of the datawide argument, and the form (wide or long) may be specified via the dataform argument.')
    if(class(try(length(dat),silent = TRUE)) == 'try-error') {
      message ('dat not specified, using datawide')
      dat <- datawide
    }
  }
  # transformedParams<-TRUE
  largeAlgebras<-TRUE
  if(nofit) fit <- FALSE
  if(fit == FALSE) carefulFit <- FALSE

  if(is.null(stationary)) stationary <- c('')
  if(all(stationary %in% 'all')) stationary<-c('T0VAR','T0MEANS','T0TIPREDEFFECT','T0TRAITEFFECT')
  
  n.latent<-ctmodelobj$n.latent
  n.manifest<-ctmodelobj$n.manifest
  Tpoints<-ctmodelobj$Tpoints
  n.TDpred<-ctmodelobj$n.TDpred
  n.TIpred<-ctmodelobj$n.TIpred
  startValues<-ctmodelobj$startValues
  
  manifestNames<-ctmodelobj$manifestNames
  latentNames<-ctmodelobj$latentNames
  TDpredNames<-ctmodelobj$TDpredNames
  TIpredNames<-ctmodelobj$TIpredNames
  
      
  
  if(dataform != 'wide' & dataform !='long') stop('dataform must be either "wide" or "long"!')
  if(dataform == 'long'){
    idcol='id'
    obsTpoints=max(unlist(lapply(unique(dat[,idcol]),function(x) 
      length(dat[dat[,idcol]==x, idcol]) )))

    if(Tpoints != obsTpoints) stop(
      paste0('Tpoints in ctmodelobj = ',Tpoints,', not equal to ', obsTpoints, ', the maximum number of rows of any subject in dat'))
    
    datawide=ctLongToWide(datalong = dat,id = 'id',
      time = 'time',manifestNames = manifestNames,TDpredNames=TDpredNames,TIpredNames = TIpredNames)
    
    datawide=suppressMessages(ctIntervalise(datawide = datawide,Tpoints = Tpoints,
      n.manifest = n.manifest,manifestNames = manifestNames,
      n.TDpred=n.TDpred,n.TIpred = n.TIpred,
      TDpredNames=TDpredNames,TIpredNames = TIpredNames))
    
    if(n.TDpred > 0){
    if(any(is.na(dat[,paste0(TDpredNames)])))  message ('Missing TD predictors found - replacing NAs with zeroes')
      datawide[,paste0(TDpredNames,'_T',rep(0:(Tpoints-1),each=n.TDpred))][is.na(datawide[,paste0(TDpredNames,'_T',rep(0:(Tpoints-1),each=n.TDpred))])] <- 0
    }
  }
  if(dataform == 'wide') datawide = dat
  
  missingManifest <- is.na(match(paste0(manifestNames, "_T0"), colnames(datawide)))
  if (any(missingManifest)) {
    stop(paste("Columns for", omxQuotes(manifestNames[missingManifest]),
      "are missing from the data frame - e.g. ", paste0(
        omxQuotes(manifestNames[missingManifest]),"_T0")))
  }
  if (length(TDpredNames)) {
    missingTD <- is.na(match(paste0(TDpredNames, "_T0"), colnames(datawide)))
    if (any(missingTD)) {
      stop(paste("Columns for", omxQuotes(TDpredNames[missingTD]),
        "are missing from the data frame"))
    }
  }
  if (length(TIpredNames)) {
    missingTI <- is.na(match(paste0(TIpredNames), colnames(datawide)))
    if (any(missingTI)) {
      stop(paste("Columns for", omxQuotes(TIpredNames[missingTI]),
        "are missing from the data frame"))
    }
  }
  
  #determine objective based on number of rows
  if(objective=='auto' & nrow(datawide) == 1) objective<-'Kalman'
  if(objective=='auto' & nrow(datawide) > 1 & n.TIpred + n.TDpred ==0 ) objective<-'mxRAM'
  if(objective=='auto' & nrow(datawide) > 1 & n.TIpred + n.TDpred >0 ) objective<-'mxRAM'
  
  #capture arguments
  ctfitargs<- list(stationary, optimizer, verbose,
    retryattempts, carefulFit, showInits, meanIntervals, fit, objective, discreteTime, asymptotes, transformedParams)
  names(ctfitargs)<-c('stationary', 'optimizer', 'verbose',
    'retryattempts', 'carefulFit', 'showInits', 'meanIntervals', 'fit', 'objective', 'discreteTime', 'asymptotes', 'transformedParams')
  
  #ensure data is a matrix
  datawide<-as.matrix(datawide)
  
  ####check data contains correct number of columns (ignore if discreteTime specified)
  neededColumns<-Tpoints * n.manifest+(Tpoints - 1)+n.TDpred * (Tpoints)+n.TIpred
  if(ncol(datawide)!=neededColumns & discreteTime != TRUE) stop("Number of columns in data (", paste0(ncol(datawide)), ") do not match model (", paste0(neededColumns), ")")
  
  if(discreteTime == TRUE){
    if(asymptotes==TRUE) stop ('Cannot estimate over asymptotic parameters for discrete time models yet!')
    message('discreteTime==TRUE set -- timing information ignored. Parameter estimates will *not* correspond with those from continous time models.')
    #     carefulFit<-FALSE
    #     stationary<-NULL    
    # if(transformedParams==TRUE){
    #   stop('For discreteTime=TRUE you must also set transformedParams=FALSE')
    # }
  }
  
  
  
  
  if(n.TDpred>0 & objective != 'Kalman' & objective != 'Kalmanmx'){ 
    #### if all tdpreds non missing, use observed covariance and means
    if(all(!is.na(datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-1), each=n.TDpred))]))){
      temp<-cov(datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-1), each=n.TDpred)),drop=FALSE]) + 
        diag(.000001, n.TDpred*(Tpoints))
      ctmodelobj$TDPREDVAR= t(chol(Matrix::nearPD(
        temp)$mat))
      ctmodelobj$TDPREDMEANS[,]=apply(datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-1), each=n.TDpred))],2,mean)
      message('No missing time dependent predictors - TDPREDVAR and TDPREDMEANS fixed to observed moments for speed')
    }

    #check for 0 variance predictors for random predictors implementation (not needed for Kalman because fixed predictors)
    ####0 variance predictor fix
    varCheck<-try(any(diag(stats::cov(datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-1), each=n.TDpred))],
      use="pairwise.complete.obs"))==0))
    if(class(varCheck)=='try-error' || any(is.na(varCheck))) {
      warning('unable to compute covariance matrix for time dependent predictors - unstable estimates may result if any variances are 0')
      varCheck<-FALSE
    }
    
    if(varCheck==TRUE &
        all(is.na(suppressWarnings(as.numeric(diag(ctmodelobj$TDPREDVAR))))) ) {
      ctmodelobj$TDPREDVAR <- diag(.1,n.TDpred*(Tpoints))
      message(paste0('Time dependent predictors with 0 variance and free TDPREDVAR matrix detected - fixing TDPREDVAR matrix diagonal to 0.01 to allow estimation.'))
    }
    
  }
  
  
  
  
  
  
  ### check single subject model adequately constrained and warn
  
  if(nrow(datawide)==1 & 'T0VAR' %in% stationary ==FALSE & 'T0MEANS' %in% stationary == FALSE & 
      all(is.na(suppressWarnings(as.numeric(ctmodelobj$T0VAR[lower.tri(ctmodelobj$T0VAR,diag=TRUE)])))) & 
      all(is.na(suppressWarnings(as.numeric(ctmodelobj$T0MEANS)))) & fit == TRUE) stop('Cannot estimate model for single individuals unless either 
        T0VAR or T0MEANS matrices are fixed, or set to stationary')
  
  if(objective == "Kalman" | objective == "Kalmanmx"){
    message('Single subject dataset or Kalman objective specified - ignoring any specified between subject 
      variance matrices (TRAITVAR, T0TRAITEFFECT,MANIFESTTRAITVAR, TIPREDEFFECT, TIPREDVAR, TDTIPREDCOV, T0TIPREDEFFECT)')
    n.TIpred<-0
  }
  
  
  if(objective=='cov') {
    if(nrow(datawide)==1) stop('Covariance based estimation impossible for single subject data')
    message('Covariance based estimation requested - meanIntervals set to TRUE')
    meanIntervals<-TRUE
  }
  
  
  
  
  ###check which extensions (e.g. traits) need to be included 
  if(any(ctmodelobj$TRAITVAR!=0) & nrow(datawide) > 1 )   traitExtension <- TRUE
  if(all(ctmodelobj$TRAITVAR==0) | nrow(datawide)==1 | objective=='Kalman' | objective=='Kalmanmx')    traitExtension <- FALSE
  
  
  if(any(ctmodelobj$MANIFESTTRAITVAR != 0)) manifestTraitvarExtension <- TRUE
  if(all(ctmodelobj$MANIFESTTRAITVAR == 0)  | nrow(datawide)==1 | objective=='Kalman'| objective=='Kalmanmx') manifestTraitvarExtension <- FALSE
  
  ## if Kalman objective, rearrange data to long format and set Tpoints to 2 (so only single discrete algebras are generated)
  if(objective=='Kalman' | objective=='Kalmanmx') {
    
    if(n.TDpred >0){
      if(any(is.na(datawide[, paste0(TDpredNames, '_T', 0:(Tpoints-1))] ))) stop('NA predictors are not possible with Kalman objective')
    }
    if(n.TIpred >0) message('Time independent predictors are not possible with single subject data, ignoring')  
    
    n.subjects<-nrow(datawide) 
    datawide<-ctWideToLong(datawide, #if using Kalman, convert data to long format
      manifestNames=manifestNames, TDpredNames=TDpredNames, TIpredNames=TIpredNames, 
      n.manifest=n.manifest, 
      Tpoints=Tpoints, 
      n.TIpred=n.TIpred, 
      n.TDpred=n.TDpred)
    
    
    colnames(datawide)[which(colnames(datawide)=='dT')]<-'dT1'
    
    if(objective == 'Kalmanmx') {
      datawide<-ctDeintervalise(datawide,dT='dT1')
      colnames(datawide)[which(colnames(datawide)=='time')] <-'dT1'
    }
    
    Tpoints<-2
    firstObsDummy<-matrix(c(1,rep(NA,times=nrow(datawide)-1)), nrow=nrow(datawide))
    for(i in 2:nrow(datawide)){
      firstObsDummy[i]<-ifelse(datawide[i,'id'] == datawide[i-1,'id'], 0, 1) #if new subject set to 1, else 0
    }
    colnames(firstObsDummy)<-'firstObsDummy'
    datawide<-cbind(datawide,firstObsDummy) 
    
    #     datawide<-rbind(c(1,rep(NA,n.manifest+n.TIpred+n.TDpred),0,1),datawide) #add empty first row so first time point included
    
  }
  
  
  #function to process ctModel specification:  seperate labels and values, fixed and free, and generate start values
  processInputMatrix <- function(x, symmetric = FALSE, diagadd = 0, randomscale=0.01, addvalues=FALSE, chol=FALSE){
    
    inputm<-x[[1]]
    
    free<-suppressWarnings(is.na(matrix(as.numeric(inputm), nrow = nrow(inputm), ncol = ncol(inputm))))
    
    labels <- ctLabel(TDpredNames=TDpredNames, TIpredNames=TIpredNames, manifestNames=manifestNames, latentNames=latentNames, matrixname=names(x), n.latent=n.latent, 
      n.manifest=n.manifest, n.TDpred=n.TDpred, n.TIpred=n.TIpred, Tpoints=Tpoints)
    
    labels[free==TRUE]<-inputm[free==TRUE]
    labels[free==FALSE]<-NA
    
    values <- matrix(round(stats::rnorm(length(inputm), 0, randomscale), 3), nrow = nrow(inputm), ncol = ncol(inputm)) #generate start values according to randomscale
    if(diagadd!= 0)  values <- values+diag(diagadd, ncol(inputm)) #increase diagonals if specified
    if(all(addvalues!=FALSE)) values <- values+addvalues #add values if specified
    if(symmetric == TRUE)  {values[row(values) > col(values)] <- t(values)[col(values) < row(values)]} #set symmetric if necessary
    values[free==FALSE] <- as.numeric(inputm[free==FALSE])
    
    if(!is.null(startValues)){ #if there are some startValues specified, check each part of x
      for( i in 1:length(startValues)){
        values[which(inputm %in% names(startValues[i]))] <- startValues[i]
      }
    }
    
    if(chol==TRUE){
      # if(any(diag(values)=='FFF0')) warning('Diagonal element of log variance input matrix fixed to 0 - to fix variance to 0, diagonal elements should tend towards -Inf, e.g. -9999',immediate. = TRUE)
      labels[row(diag(nrow(inputm))) < col(diag(nrow(inputm)))] <- NA
      values[row(diag(nrow(inputm))) < col(diag(nrow(inputm)))] <- 0
      free[row(diag(nrow(inputm))) < col(diag(nrow(inputm)))] <- FALSE
      #     labels[upper.tri(labels)]<-t(labels)[upper.tri(labels)] ### use this to generate symmetric matrices from cholesky
      #     values[upper.tri(values)]<-t(values)[upper.tri(values)]
      #     free[upper.tri(free)]<-t(free)[upper.tri(free)]
    }
    
    
    output<-list(values, labels, free)
    names(output)<-c('values', 'labels', 'free')
    return(output)
  }
  
  T0VAR <- processInputMatrix(ctmodelobj['T0VAR'], symmetric = FALSE, randomscale=0, diagadd = 10, chol=TRUE)
  
  T0MEANS <- processInputMatrix(ctmodelobj["T0MEANS"], symmetric = FALSE, randomscale=1, diagadd = 0)
  MANIFESTMEANS <- processInputMatrix(ctmodelobj["MANIFESTMEANS"], symmetric = FALSE, randomscale=1, diagadd = 0)
  LAMBDA <- processInputMatrix(ctmodelobj["LAMBDA"], symmetric = FALSE, randomscale=.1, addvalues=1, diagadd = 0)
  MANIFESTVAR <- processInputMatrix(ctmodelobj["MANIFESTVAR"],  symmetric = FALSE, randomscale=0, diagadd = 3)    
  
  DRIFT <- processInputMatrix(ctmodelobj["DRIFT"],  symmetric = FALSE,randomscale=0, 
    addvalues= ifelse(crossEffectNegStarts==TRUE,-.05,0), diagadd=ifelse(discreteTime==TRUE,.5,-.4))
  
  DIFFUSION <- processInputMatrix(ctmodelobj["DIFFUSION"], symmetric = FALSE, randomscale=0, diagadd = 10)      
  
  CINT <- processInputMatrix(ctmodelobj["CINT"], randomscale=.1)    
  
  if(transformedParams==TRUE){
    diag(T0VAR$values) <- log(diag(T0VAR$values))
    diag(T0VAR$values)[diag(T0VAR$values)== -Inf] <- -999
    
    diag(MANIFESTVAR$values) <- log(diag(MANIFESTVAR$values))
    diag(MANIFESTVAR$values)[diag(MANIFESTVAR$values)== -Inf] <- -999
    
    #     if(any(diag(DRIFT$values) >=0)) {
    #       message('transformedParams=TRUE and non negative DRIFT diagonal specified, setting to -.00001.')
    #       DRIFT$values[diag(DRIFT$values) >=0]<- -.00001
    #     }
    #     diag(DRIFT$values) <- suppressWarnings(log(-diag(DRIFT$values)) )
    # diag(DRIFT$values)[is.nan(diag(DRIFT$values)) | diag(DRIFT$values) == -Inf ] <- -999
    
    if(any(diag(DIFFUSION$values) <=0)) message('transformedParams=TRUE and non positive DIFFUSION diagonal specified, setting to .00001.')
    diag(DIFFUSION$values) <- suppressWarnings(log(diag(DIFFUSION$values)- .00001) )
    diag(DIFFUSION$values)[is.nan(diag(DIFFUSION$values)) | diag(DIFFUSION$values)== -Inf] <- -999
  }
  
  
  if(traitExtension == TRUE){ #if needed, process and include traits in matrices
    TRAITVAR <- processInputMatrix(ctmodelobj["TRAITVAR"],symmetric = FALSE, diagadd = 3, randomscale=0,chol=TRUE)
    T0TRAITEFFECT <- processInputMatrix(ctmodelobj["T0TRAITEFFECT"],symmetric = FALSE, diagadd = 3, randomscale=0,chol=FALSE)
    if(transformedParams==TRUE){
      diag(TRAITVAR$values) <- log(diag(TRAITVAR$values))
      diag(TRAITVAR$values)[diag(TRAITVAR$values)== -Inf] <- -999
    }
  }
  
  
  if(manifestTraitvarExtension == TRUE){
    MANIFESTTRAITVAR <- processInputMatrix(ctmodelobj["MANIFESTTRAITVAR"],  symmetric = FALSE, randomscale=0, diagadd = 3,chol=TRUE)
    if(transformedParams==TRUE){
      diag(MANIFESTTRAITVAR$values) <- log(diag(MANIFESTTRAITVAR$values))
      diag(MANIFESTTRAITVAR$values)[diag(MANIFESTTRAITVAR$values)== -Inf] <- -999
    }
  }
  
  
  if(traitExtension == TRUE && !is.null(ctmodelobj$TRAITTDPREDCOV)) {
    TRAITTDPREDCOV <- processInputMatrix(ctmodelobj["TRAITTDPREDCOV"], symmetric = FALSE, randomscale=0, diagadd = 0)
  }
  
  if (n.TDpred > 0){
    TDPREDMEANS <- processInputMatrix(ctmodelobj["TDPREDMEANS"], symmetric = FALSE, diagadd = 0)
    TDPREDEFFECT <- processInputMatrix(ctmodelobj["TDPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=0)
    T0TDPREDCOV <- processInputMatrix(ctmodelobj["T0TDPREDCOV"], symmetric = FALSE, diagadd = 0, randomscale=0.001) 
    
    TDPREDVAR <- processInputMatrix(ctmodelobj["TDPREDVAR"], symmetric = FALSE, diagadd = 3, randomscale=0.01,chol=TRUE) 
    if(transformedParams==TRUE){
      diag(TDPREDVAR$values) <- log(diag(TDPREDVAR$values))
      diag(TDPREDVAR$values)[diag(TDPREDVAR$values)== -Inf] <- -999
    }
  }
  
  if (n.TIpred > 0){ 
    TIPREDMEANS <- processInputMatrix(ctmodelobj["TIPREDMEANS"], symmetric = FALSE, diagadd = 0)
    TIPREDEFFECT <- processInputMatrix(ctmodelobj["TIPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=.01)
    T0TIPREDEFFECT <- processInputMatrix(ctmodelobj["T0TIPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=.01)
    TIPREDVAR <- processInputMatrix(ctmodelobj["TIPREDVAR"], symmetric = FALSE, diagadd = 3, randomscale=0,chol=TRUE)    
    if(transformedParams==TRUE){
      diag(TIPREDVAR$values) <- log(diag(TIPREDVAR$values))
      diag(TIPREDVAR$values)[diag(TIPREDVAR$values)== -Inf] <- -999
    }
  }
  
  if(n.TIpred > 0 & n.TDpred > 0) TDTIPREDCOV <- processInputMatrix(ctmodelobj["TDTIPREDCOV"], symmetric = FALSE, randomscale=0)    
  
  #     
  #     returnAllLocals()
  #### end continuous matrix section
  
  
  
  
  
  
  ###section to define base RAM matrices
  if(objective!='Kalman' & objective!='Kalmanmx') { #configure matrices 
    #   defineRAM <- function(){  
    
    #basic indexes to reuse, and modify if changed
    latentstart <- 1
    latentend <- n.latent * Tpoints
    latentExtent<-latentend #includes *all* latents - traits and manifest traits
    traitend <- latentend #because no traits yet
    latenttraitend <- traitend
    manifeststart <- n.latent * Tpoints+1
    manifestend <- manifeststart-1 + n.manifest * Tpoints
    manifestExtent <-manifestend # includes *all* manifests - predictors also
    intervalstart <- manifestend+1
    intervalend <- manifestend+Tpoints - 1
    
    
    
    #create A and S matrices   
    A<-list()
    S<-list()
    A$values <- matrix(0, nrow = manifestend, ncol = manifestend)
    A$labels <- matrix(NA, nrow = manifestend, ncol = manifestend)
    S$values <- matrix(0, nrow = manifestend, ncol = manifestend)
    S$labels <- matrix(NA, nrow = manifestend, ncol = manifestend)
    
    
    
    
    #DRIFT constraints
    diag(.8,n.latent) -> # discrete drift values goes into A$values 
      A$values[(row(A$values) - 1 - n.latent)%/%n.latent ==  # when the rows and columns, grouped by n.latent, 
          (col(A$values) - 1)%/%n.latent & row(A$values) <= latentend] #are equal, and not greater than the total number of latent variables
    
    #create expd(discrete drift) labels
    expdLabels <- cbind(paste0("discreteDRIFT_T", 
      rep(1:(Tpoints - 1), each = n.latent^2), 
      "[", seq(1, n.latent, 1), 
      ",", 
      rep(1:n.latent, each = n.latent), "]"))
    
    
    #insert discreteDRIFT_T's into A$labels    
    expdLabels-> A$labels  [(row(A$labels) - 1 - n.latent)%/%n.latent == #this groups rows by multiples of n.latent,offset by n.latent
        (col(A$labels) - 1)%/%n.latent & #this groups columns by multiples of n.latent.  when row and column groups are equal,
        row(A$labels) <= latentend] #and less than n.latent * Tpoints, then the groups are replaced by expdlabels
    
    
    
    #LAMBDA matrix
    LAMBDA$values ->  A$values[(row(A$values) - 1 - latentend)%/%n.manifest == #insert LAMBDA loadings into A$values
        (col(A$values) - 1)%/%n.latent & col(A$values) <= latentend]
    
    LAMBDA$ref <- matrix(paste0('LAMBDA[', rep(1:n.manifest, times=n.latent), ',', rep(1:n.latent, each=n.manifest), ']'), nrow=n.manifest)    
    LAMBDA$ref ->  A$labels[(row(A$labels) - 1 - latentend) %/% n.manifest == #this groups rows by multiples of n.manifest, offset by n.latent * Tpoints
        (col(A$labels) - 1)%/%n.latent & #this groups columns by multiples of n.latent.  when row and column groups are equal, 
        row(A$labels) <= manifestend] #and less than the total number of variables, then the groups are replaced by
    
    
    #measurement residuals
    diag(.8,n.manifest) ->  S$values[(row(S$values) - latentend - 1)%/%n.manifest ==  
        (col(S$values) - latentend - 1) %/% n.manifest &
        col(S$values) >  latentend &
        row(S$values) >  latentend]
    
    MANIFESTVAR$ref <- matrix(paste0('MANIFESTVAR[', 
      rep(1:n.manifest, times=n.manifest), 
      ',', 
      rep(1:n.manifest, each=n.manifest), 
      ']'), 
      nrow=n.manifest) 
    
    MANIFESTVAR$ref ->   S$labels[(row(S$labels) - latentend - 1) %/% n.manifest ==  
        (col(S$labels) - latentend - 1) %/% n.manifest &
        col(S$labels) >  latentend]
    
    
    #latent (dynamic) residuals    
    diag(1,n.latent) ->  S$values[(row(S$values) - 1) %/%n.latent == #initial DIFFUSION values with fixed coding
        (col(S$values) - 1) %/% n.latent &
        col(S$values) <= latentend]
    
    discreteDIFFUSIONlabels <- paste0("discreteDIFFUSION_T", rep(1:Tpoints - 1, each = n.latent^2), 
      "[", 
      1:n.latent, 
      ",", rep(1:n.latent,each=n.latent),"]") 
    
    
    discreteDIFFUSIONlabels-> S$labels[(row(S$labels) - 1)%/%n.latent ==  # discreteDIFFUSIONlabels goes into S$labels where the row and column groups
        (col(S$labels) - 1) %/% n.latent & # which are based on number of manifest variables, are equal, 
        col(S$labels) <= latentend] # and less than total number of latent.
    
    
    #initial latent residuals
    diag(10,n.latent) ->  S$values[1:n.latent, 1:n.latent] 
    
    T0VAR$ref <- paste0("T0VAR[", 1:n.latent, ",", rep(1:n.latent, each=n.latent), "]")
    
    T0VAR$ref -> S$labels[1:n.latent, 1:n.latent] #input initial variance refs to Slabel matrix
    
    
    ### 3. Filter matrix
    FILTER<-list()
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = n.latent * Tpoints), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(paste0(rep(latentNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.latent)), 
      paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest)))
    
    
    FILTERnamesx     <- paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest))
    
    ### 4. M matrix
    M<-list()
    M$values    <- matrix(c(T0MEANS$values, #insert appropriate T0MEANS values
      rep(T0MEANS$values, each = Tpoints-1), #follow with fixed params for later time points
      rep(MANIFESTMEANS$values, times = Tpoints)), #then insert manifest intercepts
      nrow = manifestend, ncol = 1)
    
    
    latentMlabels <- matrix(paste0("discreteCINT_T", #this name is referenced below to check which parameters to free
      0:(latentend - 1)%/%n.latent, 
      "[", 1:n.latent, ",1]"), ncol = 1)#construct continuous latent mean labels
    
    
    latentMlabels[1:n.latent, ]  <-  paste0("T0MEANS[", 1:n.latent, ",1]") #refer to T0MEANS matrix
    
    
    M$labels <- matrix(c(latentMlabels, paste0('MANIFESTMEANS',
      if('MANIFESTMEANS' %in% ctmodelobj$timeVarying) paste0('_T',rep(0:(Tpoints-1),each=n.manifest)),
      '[',rep(1:n.manifest, Tpoints),',1]')), nrow = manifestend) #insert manifest mean labels
    
    #     returnAllLocals() #return objects from this base matrices function to parent
  }#close base matrices definition function
  
  #### end RAM matrix section
  
  
  
  
  
  ###section to append RAM matrices with process traits
  if(objective!='Kalman' & traitExtension == TRUE & objective!='Kalmanmx'){ #if needed, process and include traits in matrices
    #   traitMatrices <- function(){
    
    #update indices
    manifeststart <- manifeststart+n.latent
    manifestend <- manifestend+n.latent
    traitstart <- latentend+1
    traitend <- latentend+n.latent
    latenttraitend <- traitend
    latentExtent <- traitend
    manifestExtent <- manifestend
    intervalstart <- manifestend+1
    intervalend <- manifestend+Tpoints - 1
    
    #function to insert n.latent rows and columns of single specified value into matrices 
    insertProcessTraitsToMatrix <- function(x, value){      
      x <- rbind(x[latentstart:latentend, ], #insert rows and columns for traits
        matrix(value, nrow = n.latent, ncol = ncol(x)), 
        x[(latentend+1):nrow(x), ])
      x <- cbind(x[, latentstart:latentend], 
        matrix(value, ncol = n.latent, nrow = nrow(x)), 
        x[, (latentend+1):ncol(x)])
      return(x)
    }
    
    #insert extra rows and columns to RAM matrices
    A$values <- insertProcessTraitsToMatrix(A$values, 0) 
    A$labels <- insertProcessTraitsToMatrix(A$labels, NA)  
    S$values <- insertProcessTraitsToMatrix(S$values, 0)
    S$labels <- insertProcessTraitsToMatrix(S$labels, NA)
    
    #trait loadings
    
    # #new trait loadings to manifest based on lambda
    # A$labels[manifeststart:manifestend, 
    #   (latentend+1):(latentend+n.latent)] <- paste0('LAMBDA[',
    #     rep(1:n.manifest,times=Tpoints*n.latent),',',
    #     rep(1:n.latent,each=Tpoints*n.manifest),']')
    # 
    # A$values[manifeststart:manifestend, 
    #   (latentend+1):(latentend+n.latent)] <- LAMBDA$values[cbind(
    #     rep(1:n.manifest,times=Tpoints*n.latent),
    #     rep(1:n.latent,each=Tpoints*n.manifest))]
    
    #trait loadings to latent
    # 
    if(!discreteTime) A$labels[cbind(rep(1:latentend,each=n.latent), (latentend+1):(latentend+n.latent))] <- paste0('discreteTRAIT_T',
      rep(0:(Tpoints-1),each=n.latent^2),'[',rep(1:n.latent,each=n.latent),',',1:n.latent,']')
    
    if(discreteTime) {
      for(i in 1:(Tpoints-1)){
        # 
        A$values[(i*n.latent+1):(i*n.latent+n.latent), (latentend+1):(latentend+n.latent)] <- diag(1,n.latent)
      }
    }
    
    
    #trait variance
    S$values[(latentend+1):(latentend+n.latent), (latentend+1):(latentend+n.latent)] <- diag(1,n.latent)
    TRAITVAR$ref <- matrix(paste0("TRAITVAR[", indexMatrix(symmetrical = TRUE, dimension = n.latent, sep = ","), "]"), nrow = n.latent)
    S$labels[(latentend+1):(latentend+n.latent), (latentend+1):(latentend+n.latent)] <- TRAITVAR$ref    
    
    #T0 trait effect
    A$values[1:n.latent, (latentend+1):(latentend+n.latent)] <- diag(1,n.latent)
    T0TRAITEFFECT$ref <- matrix(paste0("T0TRAITEFFECT[", indexMatrix(symmetrical = FALSE, dimension = n.latent, sep = ","), "]"), nrow = n.latent)
    A$labels[1:n.latent, (latentend+1):(latentend+n.latent)] <- T0TRAITEFFECT$ref 
    
    #M matrices
    M$values <- matrix(c(M$values[1:(n.latent * Tpoints), ], 
      rep(0, each = (n.latent)), 
      M$values[(n.latent * Tpoints+1):(n.latent * Tpoints+n.manifest * Tpoints), ]), ncol = 1)
    
    M$labels <- matrix(c(M$labels[1:(n.latent * Tpoints), ], 
      rep(NA, each = (n.latent)), 
      M$labels[(n.latent * Tpoints+1):(n.latent * Tpoints+n.manifest * Tpoints), ]), ncol = 1)
    
    
    #FILTER matrices
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = n.latent * Tpoints+n.latent), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(paste0(rep(latentNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.latent)), paste0(latentNames, "Trait"), 
      paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest)))  #subtracting manifestend from latentend gets names of manifest vars from data (because index refers to matrices)
    
    #     FILTERnamesx already created in defineRAM
    
  }#close trait matrices section
  
  
  
  
  
  
  
  ###section to append matrices with manifest trait latents
  if( objective!='Kalman' & manifestTraitvarExtension == TRUE & objective!='Kalmanmx'){
    
    #update indices
    manifeststart <- manifeststart+n.manifest #adding n.manifest latent variables to matrix indices, retaining trait indices notation rather than splitting to process and manifest
    manifestend <- manifestend+n.manifest
    manifestExtent<-manifestend
    traitstart <- latentend+1
    if(traitExtension == TRUE) traitend <- traitend + n.manifest
    if(traitExtension == FALSE) traitend <- latentend + n.manifest
    manifesttraitstart <- traitend - n.manifest + 1
    intervalstart <- manifestend+1
    intervalend <- manifestend+Tpoints - 1
    
    #function to insert n.manifest rows and columns of single specified value into matrices 
    insertManifestTraitsToMatrix <- function(x, value){      
      x <- rbind(x[latentstart:(manifesttraitstart-1), ], #insert rows and columns for manifest traits
        matrix(value, nrow = n.manifest, ncol = ncol(x)), 
        x[(manifesttraitstart):nrow(x), ])
      x <- cbind(x[, latentstart:(manifesttraitstart-1)], 
        matrix(value, ncol = n.manifest, nrow = nrow(x)), 
        x[, (manifesttraitstart):ncol(x)])
      return(x)
    }
    
    #insert new columns in matrices for manifest traits
    A$labels <- insertManifestTraitsToMatrix(A$labels, NA)  
    S$values <- insertManifestTraitsToMatrix(S$values, 0)
    S$labels <- insertManifestTraitsToMatrix(S$labels, NA)
    A$values <- insertManifestTraitsToMatrix(A$values, 0) 
    
    #manifest trait variance
    MANIFESTTRAITVAR$ref<- indexMatrix(dimension=n.manifest, symmetrical=TRUE, 
      sep=',', starttext='MANIFESTTRAITVAR[', endtext=']')
    
    S$values[manifesttraitstart:traitend, manifesttraitstart:traitend] <- diag(10,n.manifest)
    S$labels[manifesttraitstart:traitend, manifesttraitstart:traitend] <- MANIFESTTRAITVAR$ref    
    #manifest trait effect on manifests
    for(i in 0:(Tpoints-1)){
      A$values[(manifeststart+(n.manifest*i)):(manifeststart+(n.manifest* (i + 1)) - 1), 
        (manifesttraitstart):(traitend)] <- diag(n.manifest)
    }
    
    #M matrices
    M$values <- matrix(c(M$values[1:(manifesttraitstart-1), ], #include M$values before MANIFESTTRAITVARs
      rep(0, each = (n.manifest)), #fix MANIFESTTRAITVAR means to 0
      M$values[(traitend+1-n.manifest):(manifestend-n.manifest), ]), ncol = 1)#include M$values after MANIFESTTRAITVARs
    
    M$labels <- matrix(c(M$labels[1:(manifesttraitstart-1), ], 
      rep(NA, each = (n.manifest)), 
      M$labels[(traitend+1-n.manifest):(manifestend-n.manifest), ]), ncol = 1)
    
    
    #F matrices
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = traitend), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(FILTERnamesy[1:(manifesttraitstart-1)], 
      paste0(manifestNames, 'Trait'), 
      FILTERnamesy[manifesttraitstart:length(FILTERnamesy)])  #subtracting manifestend from latentend gets names of manifest vars from data (because index refers to matrices)
    
    #FILTERnamesx already created
    
  }#close manifest trait matrices
  
  
  
  
  
  
  
  ####TDpred random matrix section
  if(n.TDpred > 0 & objective!='Kalman' & objective!='Kalmanmx') {
    
    
    #function to insert rows and columns of single specified value into matrices
    insertTDpredsToMatrix <- function(target, value){
      target <- rbind(target, #insert rows and columns for traits
        matrix(value, nrow = n.TDpred * (Tpoints), ncol = ncol(target)))
      target <- cbind(target, 
        matrix(value, nrow = nrow(target), ncol = n.TDpred * (Tpoints)))
      return(target)
    }
    
    #update indices
    
    predictorTDstart <- manifestend+1
    predictorTDend <- manifestend+n.TDpred * (Tpoints)
    predictorstart<-predictorTDstart
    predictorend<-predictorTDend
    
    
    #add rows and columns to matrices
    A$values <- insertTDpredsToMatrix(A$values, 0)   
    A$labels <- insertTDpredsToMatrix(A$labels, NA)  
    S$values <- insertTDpredsToMatrix(S$values, 0)
    S$labels <- insertTDpredsToMatrix(S$labels, NA)
    
    
    #create time dependent predictor effects on processes
    TDPREDEFFECT$ref <- paste0(
      "TDPREDEFFECT", 
      "[", 
      1:n.latent, 
      ",", 
      rep(1:n.TDpred, each = n.latent), 
      "]")
    
    
    A$values[cbind(rep( 1:n.latent, times=n.TDpred*Tpoints) + n.latent*rep(0:(Tpoints-1),each=n.latent*n.TDpred), 
      rep(predictorTDstart:predictorTDend, each=n.latent))] <- TDPREDEFFECT$values #insert starting values specifying fixed for algebras
    A$labels[cbind(rep( 1:n.latent, times=n.TDpred*Tpoints) + n.latent*rep(0:(Tpoints-1),each=n.latent*n.TDpred), 
      rep(predictorTDstart:predictorTDend, each=n.latent))] <- TDPREDEFFECT$ref #insert TDPREDEFFECT algebra references to A$labels    
    
    #add cov of time dependent predictors with T0
    
    T0TDPREDCOV$ref <- paste0('T0TDPREDCOV[', 1:n.latent, ',', rep( 1:(n.TDpred*(Tpoints)), each=n.latent ), ']')
    S$values[1:n.latent, predictorTDstart:predictorTDend] <- T0TDPREDCOV$values #add starting values 
    S$values[predictorTDstart:predictorTDend, 1:n.latent] <- t(T0TDPREDCOV$values) #add starting values 
    
    S$labels[1:n.latent, predictorTDstart:predictorTDend]  <-  T0TDPREDCOV$ref #insert combined labels to S matrix
    S$labels[predictorTDstart:predictorTDend, 1:n.latent]  <-  t(S$labels[1:n.latent, predictorTDstart:predictorTDend]) #insert combined labels to S matrix
    
    TDPREDVAR$ref<-paste0('TDPREDVAR[', 1:(n.TDpred*(Tpoints)), ',', rep( 1:(n.TDpred*(Tpoints)), each=n.TDpred*(Tpoints) ), ']')
    S$values[predictorTDstart:predictorTDend, predictorTDstart:predictorTDend]  <- diag(1,n.TDpred*(Tpoints))    #insert values    
    S$labels[predictorTDstart:predictorTDend, predictorTDstart:predictorTDend ] <- TDPREDVAR$ref #insert combined labels into S matrix
    
    #introduce covariance between TDpreds and traits    
    if(traitExtension == TRUE && !is.null(ctmodelobj$TRAITTDPREDCOV)){
      TRAITTDPREDCOV$ref<-paste0('TRAITTDPREDCOV[', 1:n.latent, ',', rep( 1:(n.TDpred*(Tpoints)), each=n.latent ), ']')
      S$values[traitstart:(latentend+n.latent), predictorTDstart:predictorTDend ]  <- TRAITTDPREDCOV$values #insert starting values
      S$values[predictorTDstart:predictorTDend, traitstart:(latentend+n.latent) ]  <- t(TRAITTDPREDCOV$values)#insert symmetric starting values
      
      S$labels[traitstart:(latentend+n.latent), predictorTDstart:predictorTDend ]  <- TRAITTDPREDCOV$ref #insert combined labels to s matrix      
      S$labels[predictorTDstart:predictorTDend, traitstart:(latentend+n.latent) ]  <- t(TRAITTDPREDCOV$ref) #insert symmetric labels      
    }
    
    #Means
    
    TDPREDMEANS$ref<-matrix(paste0('TDPREDMEANS[',1:(n.TDpred*(Tpoints)),',1]'),ncol=ncol(TDPREDMEANS$labels))
    M$values  <- rbind(M$values, TDPREDMEANS$values)
    M$labels  <- rbind(M$labels, TDPREDMEANS$ref)
    
    #filter matrix    
    FILTER$values    <- cbind(matrix(0, nrow = (n.manifest * Tpoints+n.TDpred*(Tpoints)), 
      ncol = manifeststart - 1), 
      diag(1,n.manifest * Tpoints+n.TDpred*(Tpoints)))
    
    FILTERnamesy <- c(FILTERnamesy, #already specified FILTERnames
      paste0(TDpredNames, "_T", rep(0:(Tpoints-1),each=n.TDpred))) #TDpred names
    
    FILTERnamesx     <- c(FILTERnamesx,
      paste0(TDpredNames, "_T", rep(0:(Tpoints-1),each=n.TDpred)))
    
    #     returnAllLocals() #return objects from this function to parent function
  }#close TD predictor matrices function
  
  
  
  
  
  
  
  ######Time independent predictors random matrix section
  if(n.TIpred > 0 & objective != 'Kalman' & objective != 'Kalmanmx') {
    
    #function to insert rows and columns of single specified value into matrices
    insertTIpredsToMatrix <- function(target, value){
      target <- rbind(target, 
        matrix(value, nrow = n.TIpred, ncol = ncol(target)))
      target <- cbind(target, 
        matrix(value, nrow = nrow(target), ncol = n.TIpred))
      return(target)
    }
    
    #update indices
    predictorTIstart <- manifestend + n.TDpred*(Tpoints) + 1
    predictorTIend <- predictorTIstart+n.TIpred-1
    predictorend<-predictorTIend
    if(n.TDpred == 0) predictorstart<-predictorTIstart
    
    #insert rows and columns to matrices
    A$values <- insertTIpredsToMatrix(A$values, 0)
    A$labels <- insertTIpredsToMatrix(A$labels, NA)  
    S$values <- insertTIpredsToMatrix(S$values, 0)
    S$labels <- insertTIpredsToMatrix(S$labels, NA)    
    
    #Effect of TIpreds on processes
    A$values[(n.latent+1):latentend, predictorTIstart:predictorTIend] <- TIPREDEFFECT$values #add rough starting values, fixed for algebras
    
    TIPREDEFFECT$ref <- paste0( #create time Tindependent predictor algebra reference labels for A matrix 
      "discreteTIPREDEFFECT_T", 
      rep(1:(Tpoints - 1), each = n.latent), 
      "[", 
      1:n.latent, 
      ",",    
      rep(1:n.TIpred, each = (Tpoints - 1) * n.latent), 
      "]")
    
    A$labels[(n.latent+1):latentend, predictorTIstart:predictorTIend] <- TIPREDEFFECT$ref #add time Tindependent predictor algebra references to A matrix
    
    #add effect of time independent predictors on t1
    T0TIPREDEFFECT$ref <- paste0( #create time Tindependent predictor algebra reference labels for A matrix 
      "T0TIPREDEFFECT", 
      "[", 
      1:n.latent, 
      ",",    
      rep(1:n.TIpred, each=n.latent), 
      "]")
    A$values[1:n.latent, predictorTIstart:predictorTIend] <- T0TIPREDEFFECT$values
    A$labels[1:n.latent, predictorTIstart:predictorTIend ]  <-  T0TIPREDEFFECT$ref 
    
    
    #add cov between TIpreds
    TIPREDVAR$ref<-paste0('TIPREDVAR[', 1:n.TIpred, ',', rep( 1:n.TIpred, each=n.TIpred ), ']')
    S$values[predictorTIstart:predictorTIend, predictorTIstart:predictorTIend]  <- diag(10,n.TIpred)
    S$labels[predictorTIstart:predictorTIend, predictorTIstart:predictorTIend ] <- TIPREDVAR$ref 
    
    #add cov between TDpreds and TIpreds
    if(n.TDpred > 0 & n.TIpred > 0){
      
      TDTIPREDCOV$ref <- paste0('TDTIPREDCOV[',rep(1:(n.TDpred*(Tpoints)),n.TIpred),',',rep(1:n.TIpred,each=n.TDpred*(Tpoints)),']')
      S$values[predictorTDstart:predictorTDend, predictorTIstart:predictorTIend]  <- TDTIPREDCOV$values        
      S$labels[predictorTDstart:predictorTDend, predictorTIstart:predictorTIend] <- TDTIPREDCOV$ref 
      
      S$values[predictorTIstart:predictorTIend, predictorTDstart:predictorTDend]  <- t(TDTIPREDCOV$values)
      S$labels[predictorTIstart:predictorTIend, predictorTDstart:predictorTDend] <- t(TDTIPREDCOV$ref)
    }
    
    #TIpred means
    
    TIPREDMEANS$ref<-matrix(paste0('TIPREDMEANS[',1:(n.TIpred),',1]'),ncol=ncol(TIPREDMEANS$labels))
    M$values  <- rbind(M$values, TIPREDMEANS$values)
    M$labels  <- rbind(M$labels, TIPREDMEANS$ref)
    
    #Filter matrix
    FILTER$values    <- cbind(matrix(0, nrow = (n.manifest * Tpoints+n.TDpred * (Tpoints)+n.TIpred), ncol = manifeststart - 1), 
      diag(1, nrow = n.manifest * Tpoints+n.TDpred * (Tpoints)+n.TIpred))
    
    FILTERnamesy <- c(FILTERnamesy, TIpredNames)      
    FILTERnamesx     <- c(FILTERnamesx, TIpredNames)
    
    #     returnAllLocals() #return objects from this function to parent function
  }#close TI predictor matrices function
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define OpenMx RAM Algebras for the continuous time drift matrix (A), intercept (INT), and error covariance (DIFFUSION)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  
  
  if(largeAlgebras==TRUE){
    
    if(meanIntervals==TRUE) datawide[,paste0('dT', 1:(Tpoints-1))] <- 
        matrix(apply(datawide[,paste0('dT', 1:(Tpoints-1)),drop=FALSE],2,mean,na.rm=TRUE), byrow=TRUE,nrow=nrow(datawide), ncol=(Tpoints-1))
    
    uniqueintervals<-c(sort(unique(c(datawide[,paste0('dT', 1:(Tpoints-1))]))))
    
    intervalsi <- matrix(apply(datawide[,paste0('dT', 1:(Tpoints-1)),drop=FALSE], 2, 
      function(x) match(x, uniqueintervals)),ncol=(Tpoints-1))
    colnames(intervalsi)<-paste0('intervalID_T',1:(Tpoints-1))
    
    if(objective != 'cov') intervalID_T<-mxMatrix(name='intervalID_T',nrow=1,ncol=(Tpoints-1), free=FALSE,
      labels=paste0('data.intervalID_T',1:(Tpoints-1)))
    
    if(objective == 'cov') intervalID_T<-mxMatrix(name='intervalID_T',nrow=1,ncol=(Tpoints-1), free=FALSE,
      values=intervalsi[1,])
    
    datawide<-cbind(datawide,intervalsi)
    
    ######## discreteDRIFT
    #discreteDRIFTallintervals
    discreteDRIFTallintervals <- list()
    for( i in 1:length(uniqueintervals)){
      if(discreteTime==FALSE) fullAlgString <- paste0("expm(DRIFT %x%", uniqueintervals[i], ")")
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DRIFT")
      
      discreteDRIFTallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFT_i", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    #discreteDRIFTbig
    partAlgString<- paste0('discreteDRIFT_i', 1:(length(uniqueintervals)),collapse=', ')
    fullAlgString <- paste0('rbind(',partAlgString,')')
    discreteDRIFTbig <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFTbig")), 
      list(theExpression = parse(text = fullAlgString)[[1]])))
    
    #discreteDRIFTtpoints
    discreteDRIFTtpoints <- list()
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE) fullAlgString <- paste0('discreteDRIFTbig[
      ((intervalID_T[1,',i,'] -1) * nlatent + 1) : (intervalID_T[1,',i,'] * nlatent),1:nlatent]')
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DRIFT")
      
      discreteDRIFTtpoints[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFT_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    nlatent<-mxMatrix(name='nlatent',nrow=1,ncol=1,free=FALSE,values=n.latent)
    
    EXPalgs<-list(nlatent, intervalID_T, discreteDRIFTtpoints, discreteDRIFTbig, discreteDRIFTallintervals)
    
    
    
    
    
    
    ######## DRIFTHATCH
    
    #     #discreteDRIFTHATCHallintervals
    #     discreteDRIFTHATCHallintervals <- list()
    #     for( i in 1:length(uniqueintervals)){
    #       fullAlgString <- paste0("omxExponential(DRIFTHATCH %x%", uniqueintervals[i], ")")
    #       
    #       discreteDRIFTHATCHallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFTHATCH_i", i)), 
    #         list(theExpression = parse(text = fullAlgString)[[1]])))
    #     }
    #     
    #     #discreteDRIFTHATCHbig
    #     partAlgString<- paste0('discreteDRIFTHATCH_i', 1:(length(uniqueintervals)),collapse=', ')
    #     fullAlgString <- paste0('rbind(',partAlgString,')')
    #     discreteDRIFTHATCHbig <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFTHATCHbig")), 
    #       list(theExpression = parse(text = fullAlgString)[[1]])))
    #     
    #     #discreteDRIFTHATCHtpoints
    #     discreteDRIFTHATCHtpoints <- list()
    #     for( i in 1:(Tpoints-1)){
    #     fullAlgString <- paste0('discreteDRIFTHATCHbig[
    #       ((intervalID_T[1,',i,']-1) * nlatent^2 + 1) : (intervalID_T[1,',i,'] * nlatent^2), 1:(nlatent^2)]')
    #       
    #       discreteDRIFTHATCHtpoints[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFTHATCH_T", i)), 
    #         list(theExpression = parse(text = fullAlgString)[[1]])))
    #     }
    #     
    #     nlatent<-mxMatrix(name='nlatent',nrow=1,ncol=1,free=FALSE,values=n.latent)
    
    # discreteDRIFTHATCHalgs<-list(discreteDRIFTHATCHtpoints, discreteDRIFTHATCHbig, discreteDRIFTHATCHallintervals)
    #     discreteDRIFTHATCHalgs<-list(DRIFTHATCH)
    #   }
    
    
    
    
    
    ######## continuous intercept
    #discreteCINTallintervals
    discreteCINTallintervals <- list()
    for( i in 1:length(uniqueintervals)){
      if(discreteTime==FALSE & asymptotes==FALSE) fullAlgString <- 
          paste0('invDRIFT %*% (discreteDRIFT_i',i, '- II) %*% CINT')
      
      if(discreteTime==FALSE & asymptotes==TRUE) fullAlgString <- paste0('(II - discreteDRIFT_i',i, ') %*% CINT')
      
      if(discreteTime==TRUE) fullAlgString <- paste0("CINT")
      
      discreteCINTallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteCINT_i", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    
    #discreteCINTbig
    partAlgString<- paste0('discreteCINT_i', 1:(length(uniqueintervals)),collapse=', ')
    fullAlgString <- paste0('rbind(',partAlgString,')')
    discreteCINTbig <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteCINTbig")), 
      list(theExpression = parse(text = fullAlgString)[[1]])))
    
    #discreteCINTtpoints
    discreteCINTtpoints <- list()
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE) fullAlgString <- paste0('discreteCINTbig[
      ((intervalID_T[1,',i,'] -1) * nlatent + 1) : (intervalID_T[1,',i,'] * nlatent),1]')
      
      if(discreteTime==TRUE) fullAlgString <- paste0("CINT")
      
      discreteCINTtpoints[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteCINT_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    INTalgs<-list(discreteCINTtpoints, discreteCINTbig, discreteCINTallintervals)
    
    
    
    
    
    
    ######## diffusion
    
    #discreteDIFFUSIONallintervals
    discreteDIFFUSIONallintervals <- list()
    for( i in 1:length(uniqueintervals)){
      if(discreteTime==FALSE & asymptotes==FALSE) fullAlgString <- 
          # paste0("(invDRIFTHATCH %*% ((discreteDRIFTHATCH_T",i,")) - invDRIFTHATCH ) %*% rvectorize(DIFFUSION)") #optimize over continuous diffusion variance
          paste0(" asymDIFFUSION  - (discreteDRIFT_i", i, " %&% asymDIFFUSION) ") 
      
      if(discreteTime==FALSE & asymptotes==TRUE) fullAlgString <- 
          paste0("DIFFUSION  - discreteDRIFT_i", i, " %&% DIFFUSION ") 
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DIFFUSION")
      
      discreteDIFFUSIONallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDIFFUSION_i", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    
    #discreteDIFFUSIONbig
    partAlgString<- paste0('discreteDIFFUSION_i', 1:(length(uniqueintervals)),collapse=', ')
    fullAlgString <- paste0('rbind(',partAlgString,')')
    discreteDIFFUSIONbig <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDIFFUSIONbig")), 
      list(theExpression = parse(text = fullAlgString)[[1]])))
    
    #discreteDIFFUSIONtpoints
    discreteDIFFUSIONtpoints <- list()
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE) fullAlgString <- paste0('discreteDIFFUSIONbig[
      ((intervalID_T[1,',i,'] -1) * nlatent + 1) : (intervalID_T[1,',i,'] * nlatent),1:nlatent]')
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DIFFUSION")
      
      discreteDIFFUSIONtpoints[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDIFFUSION_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    Qdalgs<-list(discreteDIFFUSIONtpoints, discreteDIFFUSIONbig, discreteDIFFUSIONallintervals)
    
    
    
    
    
    
    
    
    
    
    
    
  }#end large algebras
  
  
  
  
  
  
  if(largeAlgebras==FALSE){
    if(meanIntervals==TRUE) datawide[,paste0('dT', 1:(Tpoints-1))] <- 
        matrix(apply(datawide[,paste0('dT', 1:(Tpoints-1))],2,mean,na.rm=TRUE), byrow=TRUE,nrow=nrow(datawide), ncol=(Tpoints-1))
    
    ######## discreteDRIFT
    #discreteDRIFTallintervals
    discreteDRIFTallintervals <- list()
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE) fullAlgString <- paste0("omxExponential(DRIFT %x% data.dT", i, ")")
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DRIFT")
      
      discreteDRIFTallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFT_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    EXPalgs<-list(discreteDRIFTallintervals)
    
    
    
    
    
    ######## continuous intercept
    #discreteCINTallintervals
    discreteCINTallintervals <- list()
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE & asymptotes==FALSE) fullAlgString <- 
          paste0('invDRIFT %*% (discreteDRIFT_T',i, '- II) %*% CINT')
      
      if(discreteTime==FALSE & asymptotes==TRUE) fullAlgString <- paste0('(II - discreteDRIFT_T',i, ') %*% CINT')
      
      if(discreteTime==TRUE) fullAlgString <- paste0("CINT")
      
      discreteCINTallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteCINT_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    INTalgs<-list(discreteCINTallintervals)
    
    
    
    
    
    
    ######## diffusion
    
    if(discreteTime==TRUE) asymDIFFUSIONalg <- OpenMx::mxAlgebra(name='asymDIFFUSIONalg', 
      solve(II %x% II - DRIFT %x% DRIFT ) %*%  cvectorize(DIFFUSION))
    
    #discreteDIFFUSIONallintervals
    discreteDIFFUSIONallintervals <- list()
    
    
    for( i in 1:(Tpoints-1)){
      if(discreteTime==FALSE & asymptotes==FALSE) fullAlgString <- 
          # paste0("(invDRIFTHATCH %*% ((discreteDRIFTHATCH_T",i,")) - invDRIFTHATCH ) %*% rvectorize(DIFFUSION)") #optimize over continuous diffusion variance
          paste0("asymDIFFUSION  - (discreteDRIFT_T", i, " %&% asymDIFFUSION) ") 
      
      if(discreteTime==FALSE & asymptotes==TRUE) fullAlgString <- 
          paste0(" DIFFUSION  - discreteDRIFT_T", i, " %&% DIFFUSION ") 
      
      if(discreteTime==TRUE) fullAlgString <- paste0("DIFFUSION")
      
      discreteDIFFUSIONallintervals[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDIFFUSION_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    Qdalgs<-list(discreteDIFFUSIONallintervals)
  }#end small algebras
  #end base algebra definition function
  
  
  
  
  
  
  
  ####TRAIT EXTENSION ALGEBRAS
  if(traitExtension == TRUE) {
    
    traitalgs <- list()
    for(i in 1:(Tpoints - 1)){
      
      if(discreteTime==FALSE){
        #         if(asymptotes==FALSE)
        # fullAlgString <- paste0("invDRIFT %*%   (omxExponential(DRIFT %x%", defcall[i], ") - invDRIFT)")   #optimize using continuous traitvar
        # if(asymptotes==TRUE)    
        fullAlgString <- paste0("II - discreteDRIFT_T", i) #using asymptotic trait variance
      }
      
      if(discreteTime==TRUE) fullAlgString <- paste0("II")
      
      traitalgs[[i]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteTRAIT_T", i)  ),
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    #     returnAllLocals()
  }#end Trait algebra definition function
  
  
  
  #### Predictors algebra setup
  
  if (n.TIpred > 0){ #if there are fixed time independent predictors
    discreteTIPREDEFFECTalgs <- list()
    for(j in 1:(Tpoints - 1)){
      
      if(discreteTime==FALSE){        
        #         if(asymptotes==FALSE) fullAlgString <- paste0("invDRIFT %*% (omxExponential(DRIFT %x% ", defcall[i], ") - II) %*% TIPREDEFFECT")
        #         if(asymptotes==TRUE) 
        
        
        if(asymptotes==FALSE) fullAlgString <- paste0("invDRIFT %*% (discreteDRIFT_T", j, " - II) %*% TIPREDEFFECT") 
        if(asymptotes==TRUE) fullAlgString <- paste0("(II - discreteDRIFT_T", j, ") %*% TIPREDEFFECT")
      }
      if(discreteTime==TRUE) fullAlgString <- paste0("TIPREDEFFECT") 
      
      discreteTIPREDEFFECTalgs[[j]] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteTIPREDEFFECT", "_T", j)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))  
    }
  } # end predictors model section
  
  
  
  
  
  
  ## function to expand cholesky matrices for untransformed estimation
  dechol<-function(matrixname){
    out<-list()
    for(x in c('labels','values','free')) {
      out[[x]]<-get(matrixname,envir=sys.frame(-1))[[x]]
      if(x=='values') {
        out[[x]] <- out[[x]] %*% t(out[[x]]) 
      } else 
      { out[[x]][upper.tri(out[[x]])] <- 
        out[[x]][lower.tri(out[[x]])]
      }
    }
    return(out)
  }
  
  
  
  #base model specification
  
  zerodiagnlatent <- matrix(NA, n.latent, n.latent) #set a matrix with 0's on diag for upper / lower bounds
  diag(zerodiagnlatent)<-0.000001
  zerodiagnmanifest <- matrix(NA, n.manifest, n.manifest) #set a matrix with 0's on diag for upper / lower bounds
  diag(zerodiagnmanifest)<-0.000001
  
  if('T0MEANS' %in% stationary){
    
    if(asymptotes==FALSE){
      T0MEANS$labels[T0MEANS$free==TRUE]<-paste0('asymCINT[', 1:n.latent, ',1]')[T0MEANS$free==TRUE]
      T0MEANS$free<-FALSE
      asymCINTalg<- OpenMx::mxAlgebra(name='asymCINT', -invDRIFT %*% CINT ) }
    
    if(asymptotes==TRUE) {
      T0MEANS$labels[T0MEANS$free==TRUE]<-paste0('CINT[', 1:n.latent, ',1]')[T0MEANS$free==TRUE]
      T0MEANS$free<-FALSE
    }
    
    if(discreteTime==TRUE) asymCINTalg <- OpenMx::mxAlgebra(name='asymCINT', solve(II-DRIFT) %*% CINT ) 
  }
  
  ### Set transformed or non transformed matrices
  if(discreteTime==FALSE && transformedParams==TRUE){
    DRIFT$lbound <- matrix(NA,n.latent,n.latent)
    # diag(DRIFT$lbound)<-.00001
    DRIFT$ubound<-DRIFT$lbound
    # diag(DRIFT$ubound)<- .99999
    
    DRIFT$mxmatrix<-list(
      OpenMx::mxMatrix(name = "DRIFT", type = "Full", labels = DRIFT$labels, values = DRIFT$values, free = DRIFT$free,
        ubound=DRIFT$ubound,lbound=DRIFT$lbound) 
      #       OpenMx::mxAlgebra(name='DRIFT', dimnames=list(latentNames,latentNames),
      #         negDRIFTlog - vec2diag(diag2vec(negDRIFTlog)) + vec2diag(-exp(diag2vec(negDRIFTlog))))
    )
    
    
    T0VAR$mxmatrix<-list(
      OpenMx::mxMatrix(name = "T0VARbase", values=T0VAR$values, labels=T0VAR$labels, 
        ncol=n.latent, nrow=n.latent, free=T0VAR$free, type='Full'), 
      OpenMx::mxAlgebra(name='T0VARchol', vec2diag(exp(diag2vec(T0VARbase))) + #exp of diagonal
          T0VARbase - #plus the base matrix
          vec2diag(diag2vec(T0VARbase))), #minus the diagonal of the base matrix   
      OpenMx::mxAlgebra(name='T0VAR', T0VARchol %*% t(T0VARchol))
    )
    
    
    DIFFUSION$mxmatrix<-list(
      OpenMx::mxMatrix(name = "DIFFUSIONbase", values=DIFFUSION$values, labels=DIFFUSION$labels, 
        ncol=n.latent, nrow=n.latent, free=DIFFUSION$free, type='Full'), 
      OpenMx::mxAlgebra(name='DIFFUSIONchol', vec2diag(exp(diag2vec(DIFFUSIONbase))) + #inverse log link for diagonal
          DIFFUSIONbase - #plus the base matrix
          vec2diag(diag2vec(DIFFUSIONbase))), #minus the diagonal of the base matrix   
      OpenMx::mxAlgebra(name='DIFFUSION', DIFFUSIONchol %*% t(DIFFUSIONchol))
    )
    
    MANIFESTVAR$mxmatrix<-list(
      OpenMx::mxMatrix(name = "MANIFESTVARbase", values=MANIFESTVAR$values, labels=MANIFESTVAR$labels, 
        ncol=n.manifest, nrow=n.manifest, free=MANIFESTVAR$free, type='Full'), 
      OpenMx::mxAlgebra(name='MANIFESTVARchol', vec2diag(exp(diag2vec(MANIFESTVARbase))) + #inverse log link for diagonal
          MANIFESTVARbase - #plus the base matrix
          vec2diag(diag2vec(MANIFESTVARbase))), #minus the diagonal of the base matrix   
      OpenMx::mxAlgebra(name='MANIFESTVAR', MANIFESTVARchol %*% t(MANIFESTVARchol))
    )
  }
  
  if(discreteTime==TRUE | transformedParams==FALSE){
    lboundmat=diag(1,n.latent)
    lboundmat[lboundmat==0] <- NA
    lboundmat[lboundmat==1] <- 0
    DRIFT$mxmatrix <- list( OpenMx::mxMatrix(name = "DRIFT", type = "Full", labels = DRIFT$labels, values = DRIFT$values, free = DRIFT$free))
    
    T0VAR<-dechol('T0VAR')
    DIFFUSION <- dechol('DIFFUSION')
    MANIFESTVAR <- dechol('MANIFESTVAR')
    
    T0VAR$mxmatrix<-list(
      OpenMx::mxMatrix(name = "T0VAR", values=T0VAR$values, labels=T0VAR$labels, 
        lbound=lboundmat, ncol=n.latent, nrow=n.latent, free=T0VAR$free)
    )
    
    DIFFUSION$mxmatrix<- list(
      OpenMx::mxMatrix(name = "DIFFUSION",type = "Full", labels = DIFFUSION$labels, 
        lbound=lboundmat, values = DIFFUSION$values, #DIFFUSION matrix of dynamic innovations
        free = DIFFUSION$free, nrow = n.latent, ncol = n.latent)
    )
    
    MANIFESTVAR$mxmatrix<- list(
      OpenMx::mxMatrix(name='MANIFESTVAR', free=MANIFESTVAR$free, values=MANIFESTVAR$values, 
        lbound=lboundmat, labels=MANIFESTVAR$labels, nrow=n.manifest, ncol=n.manifest)
    )
  }
  
  
  

  model  <-  OpenMx::mxModel("ctsem", #begin specifying the mxModel
    # type="RAM",
    # latentVars = paste0(rep(latentNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.latent)),
      # manifestVars = FILTERnamesx, 
      # latentVars = FILTERnamesy[!(FILTERnamesy %in% FILTERnamesx)],
    mxData(observed = datawide, type = "raw"), 
    
    mxMatrix(type = "Iden", nrow = n.latent, ncol = n.latent, free = FALSE, name = "II"), #identity matrix
    
    mxMatrix(name='LAMBDA', free=LAMBDA$free, values=LAMBDA$values, dimnames=list(manifestNames, latentNames), 
      labels=LAMBDA$labels, nrow=n.manifest, ncol=n.latent), 
    
    
    
    mxMatrix(name = "T0MEANS", free=T0MEANS$free, labels=T0MEANS$labels, #T0MEANS matrix
      values=T0MEANS$values, nrow=n.latent, ncol=1), 
    
    DRIFT$mxmatrix, MANIFESTVAR$mxmatrix, DIFFUSION$mxmatrix, T0VAR$mxmatrix,
    
    mxMatrix(type = "Full", labels = CINT$labels, values = CINT$values, 
      free = CINT$free, nrow=nrow(CINT$values), ncol=ncol(CINT$values), name = "CINT"),  #continuous intercept matrix
    INTalgs, 	#continuous intercept discrete translation algebras
    Qdalgs, # error covariance (DIFFUSION) discrete translation algebras
    mxAlgebra(name='invDRIFT', solve(DRIFT)), 
    EXPalgs #include matrix exponential algebras to equate discrete matrix to DRIFT matrix             
  )
  
  if('MANIFESTMEANS' %in% ctmodelobj$timeVarying) {
    for(occasioni in 0:(Tpoints-1)){
      
      model<-mxModel(model, 
        mxMatrix(name=paste0('MANIFESTMEANS','_T',occasioni), free=MANIFESTMEANS$free, values=MANIFESTMEANS$values, 
          labels=paste0(MANIFESTMEANS$labels,'_T',occasioni), nrow=nrow(MANIFESTMEANS$labels), ncol=ncol(MANIFESTMEANS$labels))
      )}} else {
        model<-mxModel(model, 
          mxMatrix(name='MANIFESTMEANS', free=MANIFESTMEANS$free, values=MANIFESTMEANS$values, 
            labels=MANIFESTMEANS$labels, nrow=nrow(MANIFESTMEANS$labels), ncol=ncol(MANIFESTMEANS$labels))
        )}
  
  if(discreteTime==FALSE && asymptotes==FALSE){
    
    DRIFTHATCH <- OpenMx::mxAlgebra(DRIFT %x% II + II %x% DRIFT, name = paste0("DRIFTHATCH"))
    asymDIFFUSIONalg<- OpenMx::mxAlgebra(name='asymDIFFUSIONalg', -solve(DRIFTHATCH) %*% cvectorize(DIFFUSION))
    asymDIFFUSION <- OpenMx::mxMatrix(name='asymDIFFUSION', labels=paste0('asymDIFFUSIONalg[',1:n.latent^2,',1]'),
      values=diag(10,n.latent),nrow=n.latent,ncol=n.latent)
    model<-OpenMx::mxModel(model,DRIFTHATCH,asymDIFFUSIONalg,asymDIFFUSION)
  }
  
  if('T0VAR' %in% stationary) {
    
    if(asymptotes==FALSE){
      T0VAR$labels<-paste0('asymDIFFUSIONalg[', 1:(n.latent^2), ',1]')
      T0VAR$free<-FALSE    
      
      model<-OpenMx::mxModel(model, remove=TRUE, 'T0VAR', 'T0VARchol', 'T0VARbase')
      
      model<-OpenMx::mxModel(model, 
        mxMatrix(name = "T0VAR", values=T0VAR$values, labels=T0VAR$labels, ncol=n.latent, nrow=n.latent, free=T0VAR$free)
      )
    }    
    
    if(asymptotes==TRUE){
      T0VAR$labels<-paste0('DIFFUSION[', 1:n.latent, ',', rep(1:n.latent,each=n.latent), ']')
      T0VAR$free<-FALSE
      model<-OpenMx::mxModel(model, remove=TRUE, 'T0VAR', 'T0VARchol', 'T0VARbase')
      
      model<-OpenMx::mxModel(model, 
        mxMatrix(name = "T0VAR", values=T0VAR$values, labels=T0VAR$labels, ncol=n.latent, nrow=n.latent, free=T0VAR$free)
      )
    }
    
  }
  if('T0MEANS' %in% stationary & asymptotes!=TRUE) model<-OpenMx::mxModel(model, asymCINTalg)
  
  if(objective!='Kalman' & objective != 'Kalmanmx') model<-OpenMx::mxModel(model, #include RAM matrices
    
    mxMatrix(values = A$values, free = FALSE, labels = A$labels, dimnames = list(FILTERnamesy, FILTERnamesy), name = "A"),   #directed effect matrix   
    
    mxMatrix(values = S$values, free = FALSE, labels = S$labels, dimnames = list(FILTERnamesy, FILTERnamesy), name = "S"),   #symmetric effect matrix
    
    mxMatrix(values = FILTER$values, free = FALSE, dimnames = list(FILTERnamesx, FILTERnamesy), name = "F"),  #filter matrix
    
    mxMatrix(free = FALSE, values = t(M$values), labels = t(M$labels), dimnames = list(1, FILTERnamesy), name = "M") #mean matrix
  )
  
  
  ###Trait variance 
  if(traitExtension==TRUE){
    
    if('T0TRAITEFFECT' %in% stationary){
      
      # if(asymptotes==FALSE){
      #   T0TRAITEFFECT$labels[T0TRAITEFFECT$free==TRUE] <-
      #     paste0('T0TRAITEFFECTalg[', 1:n.latent, ',', rep(1:n.latent, each=n.latent), ']')[T0TRAITEFFECT$free==TRUE]
      #   
      #   T0TRAITEFFECT$free <-FALSE
      #   T0TRAITEFFECTalg<- OpenMx::mxAlgebra(name='T0TRAITEFFECTalg', -invDRIFT)
      #   if(discreteTime==TRUE) T0TRAITEFFECTalg<- OpenMx::mxAlgebra(name='T0TRAITEFFECTalg', solve(II - DRIFT))
      # }
      # 
      # if(asymptotes==TRUE){
      T0TRAITEFFECT$labels[T0TRAITEFFECT$free==TRUE] <- NA
      T0TRAITEFFECT$values[T0TRAITEFFECT$free==TRUE] <- diag(n.latent)[T0TRAITEFFECT$free==TRUE]
      T0TRAITEFFECT$free <-FALSE
      # }
    }
    
    if(discreteTime==FALSE && transformedParams==TRUE){
      model <- OpenMx::mxModel(model, 
        traitalgs,
        OpenMx::mxMatrix(name = "TRAITVARbase", values=TRAITVAR$values, labels=TRAITVAR$labels, 
          ncol=n.latent, nrow=n.latent, free=TRAITVAR$free, type='Full'), 
        OpenMx::mxAlgebra(name='TRAITVARchol', vec2diag(exp(diag2vec(TRAITVARbase))) + #inverse log link for diagonal
            TRAITVARbase - #plus the base matrix
            vec2diag(diag2vec(TRAITVARbase))), #minus the diagonal of the base matrix   
        OpenMx::mxAlgebra(name='TRAITVAR', TRAITVARchol %*% t(TRAITVARchol)),
        
        OpenMx::mxMatrix(name = "T0TRAITEFFECT", values=T0TRAITEFFECT$values, labels=T0TRAITEFFECT$labels, 
          ncol=n.latent, nrow=n.latent, free=T0TRAITEFFECT$free, type='Full')
      )
    }
    if(discreteTime==FALSE && transformedParams==FALSE){
      model <- OpenMx::mxModel(model, 
        traitalgs,
        OpenMx::mxMatrix(name = "TRAITVAR", values=TRAITVAR$values, labels=TRAITVAR$labels, 
          ncol=n.latent, nrow=n.latent, free=TRAITVAR$free, type='Full'),
        OpenMx::mxMatrix(name = "T0TRAITEFFECT", values=T0TRAITEFFECT$values, labels=T0TRAITEFFECT$labels, 
          ncol=n.latent, nrow=n.latent, free=T0TRAITEFFECT$free, type='Full')
      )
    }
    
    if(discreteTime==TRUE | transformedParams==FALSE){
      TRAITVAR <- dechol('TRAITVAR')
      
      model <- OpenMx::mxModel(model, 
        mxMatrix( name = "TRAITVAR", type = "Full", labels = TRAITVAR$labels, values = TRAITVAR$values, 
          free = TRAITVAR$free),
        
        mxMatrix( name = "T0TRAITEFFECT", type = "Full", labels = T0TRAITEFFECT$labels, values = T0TRAITEFFECT$values, 
          free = T0TRAITEFFECT$free)
      )
    }
    
    # if('T0TRAITEFFECT' %in% stationary & asymptotes==FALSE) model<-OpenMx::mxModel(model, T0TRAITEFFECTalg)
    
  }
  
  ###Manifest Trait variance 
  
  if(manifestTraitvarExtension==TRUE){
    if(discreteTime==FALSE && transformedParams==TRUE){
      model <- OpenMx::mxModel(model, 
        OpenMx::mxMatrix(name = "MANIFESTTRAITVARbase", values=MANIFESTTRAITVAR$values, labels=MANIFESTTRAITVAR$labels, 
          ncol=n.manifest, nrow=n.manifest, free=MANIFESTTRAITVAR$free, type='Full'), 
        OpenMx::mxAlgebra(name='MANIFESTTRAITVARchol', vec2diag(exp(diag2vec(MANIFESTTRAITVARbase))) + #inverse log link for diagonal
            MANIFESTTRAITVARbase - #plus the base matrix
            vec2diag(diag2vec(MANIFESTTRAITVARbase))), #minus the diagonal of the base matrix   
        OpenMx::mxAlgebra(name='MANIFESTTRAITVAR', MANIFESTTRAITVARchol %*% t(MANIFESTTRAITVARchol))
      )
    }
    
    if(discreteTime==TRUE | transformedParams==FALSE){
      MANIFESTTRAITVAR <- dechol('MANIFESTTRAITVAR')
      model <- OpenMx::mxModel(model, 
        mxMatrix(type = "Full", 
          labels = MANIFESTTRAITVAR$labels, 
          values = MANIFESTTRAITVAR$values,  
          free = MANIFESTTRAITVAR$free,
          name = "MANIFESTTRAITVAR")
      )
    }
  }
  
  #model options
  originaloptimizer<- OpenMx::mxOption(NULL, "Default optimizer")
  OpenMx::mxOption(NULL, "Default optimizer", optimizer)
  
  #end base model spec
  
  
  
  
  
  
  
  
  ###predictor extension model specification
  if(n.TDpred + n.TIpred > 0){
    
    if(n.TDpred > 0 ) { #if there are fixed TD predictors 
      
      if('T0TDPREDCOV' %in% stationary){
        T0TDPREDCOV$values[T0TDPREDCOV$free==TRUE] <-0
        T0TDPREDCOV$free<-FALSE     
      }
      
      model <- OpenMx::mxModel(model, 
        mxMatrix(type = "Full", 
          labels = TDPREDEFFECT$labels, values = TDPREDEFFECT$values, free = TDPREDEFFECT$free, name = "TDPREDEFFECT"))
      
      if(objective!='Kalman' & objective != 'Kalmanmx'){
        
        if(discreteTime==FALSE && transformedParams==TRUE){
          model <- OpenMx::mxModel(model, 
            OpenMx::mxMatrix(name = "TDPREDVARbase", values=TDPREDVAR$values, labels=TDPREDVAR$labels, 
              ncol=n.TDpred*(Tpoints), nrow=n.TDpred*(Tpoints), free=TDPREDVAR$free, type='Full'), 
            OpenMx::mxAlgebra(name='TDPREDVARchol', vec2diag(exp(diag2vec(TDPREDVARbase))) + #inverse log link for diagonal
                TDPREDVARbase - #plus the base matrix
                vec2diag(diag2vec(TDPREDVARbase))), #minus the diagonal of the base matrix   
            OpenMx::mxAlgebra(name='TDPREDVAR', TDPREDVARchol %*% t(TDPREDVARchol))
          )
        }
        
        
        if(discreteTime==TRUE | transformedParams==FALSE){
          TDPREDVAR <- dechol('TDPREDVAR')
          model <- OpenMx::mxModel(model, 
            mxMatrix(type = "Full", 
              labels = TDPREDVAR$labels, values = TDPREDVAR$values, free = TDPREDVAR$free, 
              ncol=n.TDpred*(Tpoints), nrow=n.TDpred*(Tpoints),name = "TDPREDVAR")
          )
        }
        
        
        model <- OpenMx::mxModel(model, 
          mxMatrix(name='TDPREDMEANS', type='Full', labels=TDPREDMEANS$labels, free=TDPREDMEANS$free,
            values=TDPREDMEANS$values,ncol=ncol(TDPREDMEANS$labels),nrow=nrow(TDPREDMEANS$labels)),
          
          mxMatrix(type = "Full", 
            labels = T0TDPREDCOV$labels, values = T0TDPREDCOV$values, free = T0TDPREDCOV$free, name = "T0TDPREDCOV")
        )
        
        if(traitExtension==TRUE) model<-OpenMx::mxModel(model,        
          mxMatrix(labels = TRAITTDPREDCOV$labels, values = TRAITTDPREDCOV$values, free = TRAITTDPREDCOV$free, name = "TRAITTDPREDCOV")
        )
      }
    }
    
    
    if (n.TIpred > 0){ #if there are fixed time independent predictors
      
      if('T0TIPREDEFFECT' %in% stationary){
        if(asymptotes==FALSE){
          T0TIPREDEFFECT$labels[T0TIPREDEFFECT$free==TRUE] <-
            paste0('asymTIPREDEFFECT[', 1:n.latent, ',', rep(1:n.TIpred, each=n.latent), ']')[T0TIPREDEFFECT$free==TRUE]
          T0TIPREDEFFECT$free<-FALSE
          asymTIPREDEFFECTalg<- OpenMx::mxAlgebra(name='asymTIPREDEFFECT', -invDRIFT %*% TIPREDEFFECT)  
          if(discreteTime==TRUE) asymTIPREDEFFECTalg<- OpenMx::mxAlgebra(name='asymTIPREDEFFECT', (solve(II - DRIFT) %*% TIPREDEFFECT))  
          model<-OpenMx::mxModel(model, asymTIPREDEFFECTalg)
          
        }
        
        if(asymptotes==TRUE){
          T0TIPREDEFFECT$labels[T0TIPREDEFFECT$free==TRUE] <-
            paste0('TIPREDEFFECT[', 1:n.latent, ',', rep(1:n.TIpred, each=n.latent), ']')[T0TIPREDEFFECT$free==TRUE]
          T0TIPREDEFFECT$free<-FALSE
        }
        
      }
      
      if(discreteTime==FALSE && transformedParams==TRUE){
        model <- OpenMx::mxModel(model, 
          OpenMx::mxMatrix(name = "TIPREDVARbase", values=TIPREDVAR$values, labels=TIPREDVAR$labels, 
            ncol=n.TIpred, nrow=n.TIpred, free=TIPREDVAR$free, type='Full'), 
          OpenMx::mxAlgebra(name='TIPREDVARchol', vec2diag(exp(diag2vec(TIPREDVARbase))) + #inverse log link for diagonal
              TIPREDVARbase - #plus the base matrix
              vec2diag(diag2vec(TIPREDVARbase))), #minus the diagonal of the base matrix   
          OpenMx::mxAlgebra(name='TIPREDVAR', TIPREDVARchol %*% t(TIPREDVARchol))
        )
      }
      
      if(discreteTime==TRUE | transformedParams==FALSE){
        
        TIPREDVAR <- dechol('TIPREDVAR') 
        model <- OpenMx::mxModel(model, 
          mxMatrix(type = "Full", labels = TIPREDVAR$labels, values = TIPREDVAR$values, free = TIPREDVAR$free, name = "TIPREDVAR")
        )
      }
      
      model <- OpenMx::mxModel(model, 
        mxMatrix(type = "Full", labels = TIPREDEFFECT$labels, values = TIPREDEFFECT$values, free = TIPREDEFFECT$free, name = "TIPREDEFFECT"), 
        mxMatrix(name='T0TIPREDEFFECT', values = T0TIPREDEFFECT$values, nrow=nrow(T0TIPREDEFFECT$values), ncol=ncol(T0TIPREDEFFECT$values), 
          free=T0TIPREDEFFECT$free, labels=T0TIPREDEFFECT$labels), 
        
        mxMatrix(name='TIPREDMEANS', type='Full', labels=TIPREDMEANS$labels, free=TIPREDMEANS$free,
          values=TIPREDMEANS$values,ncol=ncol(TIPREDMEANS$labels),nrow=nrow(TIPREDMEANS$labels)),
        
        discreteTIPREDEFFECTalgs
      )
      if(n.TDpred > 0 ) model<-OpenMx::mxModel(model, 
        mxMatrix(name='TDTIPREDCOV', type='Full', labels=TDTIPREDCOV$labels, free=TDTIPREDCOV$free,
          values=TDTIPREDCOV$values,ncol=ncol(TDTIPREDCOV$labels),nrow=nrow(TDTIPREDCOV$labels))
      )
    }
    
    
    
    
    #     if(randomPredictors==FALSE & n.TDpred > 0){ #if we want to insert TD predictors via mean inputs controlled by definition variables
    #       
    #       model <- OpenMx::mxModel(model, #add TDPREDEFFECT matrix to model
    #         mxMatrix(type = "Full", 
    #           labels = TDPREDEFFECTlabels, 
    #           values = TDPREDEFFECT$values, 
    #           free = TDPREDEFFECT$free, 
    #           name = "TDPREDEFFECT"), 
    # #         T0TDPREDCOVmatrices, 
    #         TDpreddefs, TDpreddefsfull, INTalgs 
    #         ) #add TDpreddefs and intalgs and T0cov to model
    #       
    #       
    #             model <- OpenMx::mxModel(model, remove=TRUE, 'T0MEANS')        
    #             model <- OpenMx::mxModel(model, 
    #               mxAlgebra(name = "T0MEANS", T0MEANBASE + T0TDPREDCOV %*% t(TDpreddefsfull)), 
    #               mxMatrix(name='T0TDPREDCOV', nrow=nrow(T0TDPREDCOV), ncol=ncol(T0TDPREDCOV), 
    #                 free=T0TDPREDCOV$free, labels=T0TDPREDCOVlabels), 
    #               mxMatrix(name = "T0MEANBASE", free=T0MEANS$free, labels=T0MEANS$labels, 
    #                 values=T0MEANS$valuesnrow=n.latent, ncol=1)
    #             )
    #     }# end fixed TD predictors
    #     
    #     
    #     
    #     
    #     if(randomPredictors==FALSE & n.TIpred >0){ #if we have TI predictors
    #       
    #       model <- OpenMx::mxModel(model, remove=TRUE, 'T0MEANS')        
    #       model <- OpenMx::mxModel(model, 
    #         mxMatrix(type = "Full", labels = TIPREDEFFECTlabels, values = TIPREDEFFECT$values free = TIPREDEFFECT$free, name = "TIPREDEFFECT"), 
    #         INTalgs, 
    #         mxMatrix(name="TIpreddefs", labels=paste0("data.", TIpredNames), nrow=1, ncol=n.TIpred), #create a definition variable matrix
    #         mxAlgebra(name = "T0MEANS", T0MEANBASE + T0TIPREDEFFECT %*% t(TIpreddefs)), 
    #         mxMatrix(name='T0TIPREDEFFECT', nrow=nrow(T0TIPREDEFFECT), ncol=ncol(T0TIPREDEFFECT), 
    #           free=T0TIPREDEFFECT$free, labels=T0TIPREDEFFECTlabels), 
    #         mxMatrix(name = "T0MEANBASE", free=T0MEANS$free, labels=T0MEANS$labels
    #           values=T0MEANS$valuesnrow=n.latent, ncol=1)
    #       )
    #     } # end fixed TI predictor / intercept algebra
    
    #     returnAllLocals()
  } # end predictors model section
  
  
  
  
  
  
  
  #   setobjective<-function(){
  ######### objective functions
  if(objective == "mxRAM") {  
    # browser()
    model <- OpenMx::mxModel(model, 
      # type='RAM',
      # manifestVars = FILTERnamesx,
      # latentVars = FILTERnamesy[!(FILTERnamesy %in% FILTERnamesx)],
      # mxMatrix(values = FILTER$values, free = FALSE, dimnames = list(FILTERnamesx, FILTERnamesy), name = "F"),
      #       mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
      #       mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
      #       mxMatrix(type = "Iden", nrow = nrow(A$labels), ncol = ncol(A$labels), name = "bigI"), 
      mxExpectationRAM(M = "M",thresholds=ifelse(is.null(ctmodelobj$ordNames),NA,'thresh')),
      mxFitFunctionML(vector=FALSE)
    )
  }
  
  
  if(objective == "mxFIML" | objective == 'mxRowFIML') { #split RAM style matrices into components for faster processing
    
    manifestExtent<-manifestend #update manifest extents - perhaps requires predictors
    latentExtent<-latentend
    if(n.TDpred + n.TIpred > 0) manifestExtent <- manifestExtent + n.TDpred *(Tpoints-1) + n.TIpred
    if(n.TDpred + n.TIpred > 0) latentExtent <- latentExtent + n.TDpred *(Tpoints-1) + n.TIpred
    if(traitExtension==TRUE) latentExtent<-latentExtent+n.latent
    if(manifestTraitvarExtension==TRUE) latentExtent <- latentExtent + n.manifest
    
    latents<-1:(n.latent*Tpoints)
    if(traitExtension==TRUE) latents<-c(latents,(n.latent*Tpoints+1):(n.latent*Tpoints+n.latent))
    if(manifestTraitvarExtension==TRUE) latents<-c(latents,(manifesttraitstart):(traitend))
    if(n.TDpred > 0) latents<-c(latents,predictorTDstart:predictorTDend)
    if(n.TIpred >0) latents <- c(latents,predictorTIstart:predictorTIend)
    
    Amanifestvalues <- A$values[manifeststart:manifestExtent, latents]
    Smanifestvalues<-S$values[manifeststart:manifestExtent, manifeststart:manifestExtent]
    Smanifestlabels<-S$labels[manifeststart:manifestExtent, manifeststart:manifestExtent]
    Smanifestfree<-FALSE
    
    Amanifestcovvalues <- A$values[manifeststart:manifestExtent, latents] #need this to incorporate predictor and manifest covariance
    if(n.TDpred+n.TIpred > 0) {
      Amanifestcovvalues[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
        (traitend+1):(traitend+n.TIpred+n.TDpred*(Tpoints-1))] [col(diag(n.TIpred+n.TDpred*(Tpoints-1))) == 
            row(diag(n.TIpred+n.TDpred*(Tpoints-1)))] <- 1
      
      Smanifestvalues[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
        (n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1))]  <-0
      Smanifestlabels[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
        (n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1))]  <-NA
      Smanifestfree  <-FALSE
    }
    
    #base model components
    model <- OpenMx::mxModel(model, 
      
      mxMatrix(name='Alatent', values=A$values[latents, latents], 
        labels=A$labels[latents, latents], 
        free=FALSE, 
        nrow=latentExtent, ncol=latentExtent), 
      
      mxMatrix(name='Slatent', values=S$values[latents, latents], 
        labels=S$labels[latents, latents], 
        free=FALSE, 
        nrow=latentExtent, ncol=latentExtent), 
      
      mxMatrix(name='Mlatent', values=M$values[latents, 1], 
        labels=M$labels[latents, 1], 
        free=FALSE, 
        nrow=1, ncol=latentExtent), 
      
      mxMatrix(name='Smanifest', labels=Smanifestlabels, #includes predictors!!
        values=Smanifestvalues, 
        free=FALSE, 
        nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred), 
      #           nrow=n.manifest*Tpoints, ncol=n.manifest*Tpoints), 
      
      mxMatrix(name='Mmanifest', labels=M$labels[manifeststart:manifestExtent], 
        values=M$values[manifeststart:manifestExtent], free=FALSE, 
        nrow=1, ncol=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred), 
      #           nrow=1, ncol=n.manifest*Tpoints), 
      
      mxMatrix(name='Amanifest', values=Amanifestvalues, 
        labels=A$labels[manifeststart:manifestExtent, latents], 
        free=FALSE, 
        nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=latentExtent), 
      #           nrow=n.manifest*Tpoints, ncol=latentExtent), 
      
      mxMatrix(name='Amanifestcov', values=Amanifestcovvalues, 
        labels=A$labels[manifeststart:manifestExtent, latents], 
        free=FALSE, 
        nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=latentExtent), 
      
      mxMatrix(name='Ilatent', type='Iden', nrow=latentExtent, ncol=latentExtent), 
      
      mxAlgebra(name='invIminusAlatent', -solve(-Ilatent + Alatent))
    )
    
    fullSlatentlist<-'Slatent'   #begin a list of entities that sum together to generate the S matrix (ie traits, predictors, add to this)
    fullSmanifestlist<-'Smanifest'
    fullMlatentlist<-'Mlatent'
    
    fullSlatent<-eval(parse(text=paste0("mxAlgebra(name='fullSlatent', ", paste(fullSlatentlist, collapse=' + '), ")")))
    fullSmanifest<-eval(parse(text=paste0("mxAlgebra(name='fullSmanifest', ", paste(fullSmanifestlist, collapse=' + '), ")")))
    fullMlatent<-eval(parse(text=paste0("mxAlgebra(name='fullMlatent', ", paste(fullMlatentlist, collapse=' + '), ")")))
    
    model <- OpenMx::mxModel(model, 
      fullSlatent, fullSmanifest, fullMlatent, 
      mxAlgebra(Amanifestcov %&% (invIminusAlatent %&% fullSlatent) +Smanifest, name = "expCov"), 
      mxAlgebra(t(Amanifest %*% (invIminusAlatent %*% t(fullMlatent))) + Mmanifest, name = "expMean")) 
    
    #         mxMatrix(type = "Iden", nrow = nrow(A$labels, ncol = ncol(A$labels, name = "bigI"), 
    if(objective=='mxFIML'){
      model<-OpenMx::mxModel(model,
        mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx,thresholds=ifelse(is.null(ctmodelobj$ordNames),NA,'thresh')), 
        mxFitFunctionML())
    }
    
    if(objective=='mxRowFIML'){
      model<-OpenMx::mxModel(model,
        # mxMatrix(type = "Iden", nrow = n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, 
        # ncol = n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, name = "bigI"),
        
        mxAlgebra(expression=omxSelectRowsAndCols(expCov, existenceVector), name="filteredExpCov"),
        mxAlgebra(expression=omxSelectCols(expMean, existenceVector), name="filteredExpMean"),
        # mxAlgebra(expression=omxSelectRowsAndCols(bigI, existenceVector), name="filteredbigI"),
        mxAlgebra(expression= chol(filteredExpCov), name="filteredExpCovchol"),
        mxAlgebra(expression= solve(filteredExpCovchol), name="filteredExpCovcholinv"),
        mxAlgebra(expression=sum(existenceVector), name="numVar_i"),      
        mxAlgebra(expression = log(2*pi), name = "log2pi"),
        # mxAlgebra(expression=log2pi %*% numVar_i + log(det(filteredExpCov)), name ="firstHalfCalc"),
        mxAlgebra(expression=log2pi %*% numVar_i + log((det(filteredExpCovchol)^2)), name ="firstHalfCalc"),
        # mxAlgebra((filteredDataRow - filteredExpMean) %&% solve(filteredExpCov), name = "secondHalfCalc"),
        mxAlgebra((filteredDataRow - filteredExpMean) %&% (filteredExpCovcholinv %*% t(filteredExpCovcholinv)) , name = "secondHalfCalc"),
        mxAlgebra(expression=(firstHalfCalc + secondHalfCalc),name="rowAlgebra"),
        mxAlgebra(expression=sum(rowResults),name = "reduceAlgebra"),
        mxFitFunctionRow(rowAlgebra='rowAlgebra', reduceAlgebra='reduceAlgebra', 
          existenceVector='existenceVector', dimnames=FILTERnamesx)     
      )}
  }#end FIML objectives
  
  if(objective=='cov'){
    
    manifests <- c(paste0(manifestNames,'_T',rep(0:(Tpoints-1),each=n.manifest)), 
      if(n.TDpred > 0) {paste0(TDpredNames,'_T',rep(0:(Tpoints-2), times=n.TDpred))}, 
      if(n.TIpred > 0) {paste0(TIpredNames) } )
    
    covData<-matrix(Matrix::nearPD(stats::cov(datawide[,manifests],use='pairwise.complete.obs'))[["mat"]],nrow=length(manifests),
      dimnames = list(manifests,manifests))
    
    meanData <- apply(datawide[,manifests,drop=FALSE],2,mean,na.rm=TRUE)
    model<-OpenMx::mxModel(model, mxData(covData, type='cov', means=meanData, numObs=nrow(datawide)),
      mxExpectationRAM(M='M',thresholds=ifelse(is.null(ctmodelobj$ordNames),NA,'thresh')),
      mxFitFunctionML()
    )
  }
  
  if(objective=='Kalman'){ 
    
    # discreteDIFFUSIONmatrix<-OpenMx::mxMatrix(name='discreteDIFFUSIONmatrix',
    #   labels=paste0('discreteDIFFUSION_T1[',1:n.latent,',', rep(1:n.latent,each=n.latent),']'),
    #   nrow=n.latent,ncol=n.latent,free=FALSE)
    
    D<-OpenMx::mxMatrix(name='D', values=MANIFESTMEANS$values, labels=MANIFESTMEANS$labels, 
      nrow=n.manifest, ncol=1, free=MANIFESTMEANS$free)
    
    intercept <- OpenMx::mxAlgebra(name='intercept',discreteCINT_T1)
    
    u<-OpenMx::mxMatrix(name='u', values=1, nrow=1, ncol=1, free=FALSE)
    
    discreteDRIFT<-mxAlgebra(name='discreteDRIFT',  (1-firstObsDummy) %x% discreteDRIFT_T1)
    
    kalmanT0MEANS <- OpenMx::mxAlgebra(name='kalmanT0MEANS',T0MEANS)
    
    x0 <-  OpenMx::mxAlgebra(name='x0', T0MEANS)
    
    # if(n.subjects==1){ #then simple kalman
    #   
    #   model<-OpenMx::mxModel(model, 
    #     D, u,  discreteDIFFUSIONmatrix,
    #     mxExpectationStateSpace(A='discreteDRIFT_T1', B='discreteCINT_T1', C='LAMBDA', 
    #       D="D", Q='discreteDIFFUSIONmatrix', R='MANIFESTVAR', x0='T0MEANS', P0='T0VAR', u="u"), 
    #     mxFitFunctionML()
    #   )
    # }
    # 
    # if(n.subjects > 1){ #then account for further first time point observations
    
    
    #free intercepts
    #      randIntercepts<- OpenMx::mxMatrix(type = "Full", 
    #         labels = paste0('s',rep(unique(datawide[,'id']),each=n.latent),'_', CINT$labels), 
    #         values = CINT$values, 
    #         free = CINT$free, nrow=n.latent, ncol=n.subjects, name = "randIntercepts")
    #       
    #         model<-OpenMx::mxModel(model,'CINT',remove=TRUE)
    #        model<-OpenMx::mxModel(model,
    #          randIntercepts,
    #          # mxAlgebra(name='tCINTmatrix',t(CINTmatrix)),
    #          mxAlgebra(name='CINTalg',randIntercepts[,data.id]),
    #          mxMatrix(name='CINT',labels=paste0('CINTalg[',1:n.latent,',1]'),nrow=n.latent,ncol=1,type='Full')
    #        )
    # 
    # } #end multi subject kalman
    
    if(n.TDpred>0){
      interceptCINTlabels<-matrix(paste0('discreteCINT_T1[', 1:n.latent, ',','1]'), nrow=n.latent)
      interceptTDPREDlabels<-matrix(paste0('TDPREDEFFECT[', 1:n.latent, ',', rep(1:n.TDpred, each=n.latent), ']'), nrow=n.latent)
      
      model<-OpenMx::mxModel(model,
        OpenMx::mxMatrix(type='Full',name='TDPREDS',
          labels=paste0('data.', TDpredNames),free=FALSE,ncol=1,nrow=n.TDpred))
      
      intercept<-  OpenMx::mxMatrix(type='Full',name='intercept', free=FALSE , nrow=n.latent, ncol=n.TDpred+1, 
        labels=cbind(interceptCINTlabels, interceptTDPREDlabels))
      
      kalmanT0MEANS <- OpenMx::mxMatrix(type='Full',name='kalmanT0MEANS',free=FALSE,nrow=n.latent,ncol=n.TDpred+1,
        labels=c(paste0('T0MEANS[',1:n.latent,',1]'),rep(NA,n.TDpred*n.latent)),values=0)
      
      kalmanT0MEANS <- OpenMx::mxAlgebra(name='kalmanT0MEANS',cbind(T0MEANS,TDPREDEFFECT))
      
      x0 <-  OpenMx::mxAlgebra(name='x0', T0MEANS) # + TDPREDEFFECT %*% TDPREDS
      
      D <-  OpenMx::mxMatrix(name='D', nrow=n.manifest, ncol=1+n.TDpred, 
        free=c(MANIFESTMEANS$free, rep(FALSE, n.TDpred*n.manifest)), 
        values=c(MANIFESTMEANS$values, rep(0, n.TDpred*n.manifest)), 
        labels=c(MANIFESTMEANS$labels, rep(NA, n.TDpred*n.manifest)))
      
      u <-  OpenMx::mxMatrix(name='u', ncol=1, nrow=n.TDpred+1, free=FALSE, 
        values=c(1, rep(0, n.TDpred)), 
        labels=c(NA, paste0('data.', TDpredNames)))
      
      
    }
    
    
    
    model<-OpenMx::mxModel(model, 
      D, u, discreteDRIFT,intercept,kalmanT0MEANS,x0,
      
      mxAlgebra(name='B', (1-firstObsDummy) %x% intercept + firstObsDummy %x% kalmanT0MEANS), # 
      
      mxAlgebra(name='discreteDIFFUSIONwithdummy',  (1-firstObsDummy) %x% discreteDIFFUSION_T1 + firstObsDummy %x% (T0VAR)),
      
      mxMatrix(name='firstObsDummy', free=FALSE, labels='data.firstObsDummy', nrow=1, ncol=1),
      
      mxExpectationStateSpace(A='discreteDRIFT', B='B', C='LAMBDA', 
        D="D", Q='discreteDIFFUSIONwithdummy', R='MANIFESTVAR', x0='x0', P0='T0VAR', u="u",
        thresholds=ifelse(is.null(ctmodelobj$ordNames),NA,'thresh')), 
      
      mxFitFunctionML()
      
    )
    
  }
  
  
  
  
  if(objective=='Kalmanmx'){ 
    
    if(n.TDpred > 0) Bmat<-mxMatrix(name='B',values=cbind(CINT$values,TDPREDEFFECT$values),
      labels=cbind(CINT$labels,TDPREDEFFECT$labels),free=cbind(CINT$free,TDPREDEFFECT$free),
      nrow=n.latent,ncol=1+n.TDpred)
    
    if(n.TDpred == 0) Bmat<-mxMatrix(name='B',values=cbind(CINT$values),
      labels=cbind(CINT$labels),free=cbind(CINT$free),
      nrow=n.latent,ncol=1)
    
    model<-OpenMx::mxModel(model,
      
      mxMatrix(name='D', values=cbind(MANIFESTMEANS$values,matrix(0,nrow=n.manifest,ncol=n.TDpred)), 
        labels=cbind(MANIFESTMEANS$labels,matrix(NA,nrow=n.manifest,ncol=n.TDpred)),
        nrow=n.manifest, ncol=1+n.TDpred, free=cbind(MANIFESTMEANS$free,matrix(FALSE,nrow=n.manifest,ncol=n.TDpred))), 
      
      mxMatrix(name='u',values=c(1,n.TDpred),labels=c(NA,if(n.TDpred > 0) paste0('data.',TDpredNames)),nrow=1+n.TDpred,ncol=1),
      
      Bmat,
      
      mxMatrix('Full', 1, 1, name='time', labels='data.dT1'),
      
      mxExpectationSSCT(A='DRIFT', B='B', C='LAMBDA', 
        D="D", Q='DIFFUSION', R='MANIFESTVAR', x0='T0MEANS', P0='T0VAR', u="u", t='time',
        thresholds=ifelse(is.null(ctmodelobj$ordNames),NA,'thresh')), 
      
      mxFitFunctionML(vector=FALSE)
    )
  }
  
  
  
  
  
  
  #     if(objective == "mxFIMLpenalised") {  ## attempt for multigroup
  #       
  #       if(traitExtension==TRUE) penalties <- mxAlgebra(name='penalties', sum(DRIFT*DRIFT + TRAITVAR*TRAITVAR + 
  #           MANIFESTVAR*MANIFESTVAR + tr(DRIFT)/2) * FIMLpenaltyweight)
  #       if(traitExtension==FALSE) penalties <- mxAlgebra(name='penalties', sum(DRIFT*DRIFT + 
  #           MANIFESTVAR*MANIFESTVAR + tr(DRIFT)/2) * FIMLpenaltyweight)
  #       
  #       
  #       model <- OpenMx::mxModel(model, 
  #         #                    fitmodel, 
  #         mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx), 
  #         mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
  #         mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
  #         mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI"), 
  #         #                         mxData(observed = datawide, type = "raw"), 
  #         mxMatrix(type='Full', name='FIMLpenaltyweight', nrow=1, ncol=1, values=FIMLpenaltyweight, free=FALSE), 
  #         mxFitFunctionML(vector=TRUE), 
  #         penalties, 
  #         mxAlgebra(sum(fitfunction)+penalties, name='m2ll'), #
  #         mxFitFunctionAlgebra('m2ll')
  #       )
  #     }
  
  
  #       
  #       if(objective == "mxRowFIML") {
  #         model <- addmxrowobjective(model, dimlabels = FILTERnamesx)
  #         model <- OpenMx::mxModel(model, 
  #           mxRAMObjective("A", "S", "F", "M"), 
  #           mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
  #           mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
  #           mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI")
  #         )}
  #       
  #       if(objective == "CondLogLik") {
  #         model <- OpenMx::mxModel(model, type = "default", mxFitFunctionR(CondLogLik), 
  #           mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
  #           mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
  #           mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI")
  #         )}
  
  
  
  #     returnAllLocals()
  #   } #end setobjective function
  
  
  
  
  
  
  #   objective <- 'mxRAM' #only option for the moment so not included in initial arguments
  #   if(carefulFit==TRUE){
  #     targetObjective <- objective
  #     objective <- 'mxFIMLpenalised'
  #     
  #   }
  
  
  model <- OpenMx::omxAssignFirstParameters(model, indep = FALSE) #randomly selects inits for labels with 2 or more starting values
  if(showInits==TRUE & carefulFit==TRUE) { #
    message('Starting values: ')
    newstarts <- OpenMx::omxGetParameters(model)
    message(paste(names(newstarts), ": ", newstarts, "\n"))
  }
  
  #specfic stationarity constraints
  if(!('T0VAR' %in% stationary) && 'stationary' %in% T0VAR$labels ){
    model <- mxModel(model,'T0VAR',remove = TRUE)
    model <- mxModel(model,
      mxAlgebra(T0VARchol %*% t(T0VARchol),name='T0VARalg'),
      mxMatrix(name='T0VAR',free=FALSE,values=0,nrow=n.latent,ncol=n.latent,
        labels=paste0('T0VARalg[',1:n.latent,',',rep(1:n.latent,each=n.latent),']'))
    )
    
    for(i in 1:n.latent){ #covariance between all dynamic latents
      for(j in 1:n.latent){
        if( rl(T0VAR$labels[i,j]=='stationary')){
          model$T0VARbase$free[i, j] <- FALSE
          model$T0VAR$labels[i, j] <- paste0('asymDIFFUSION[',i,',',j,']')
        }
      }
    }
  }
  
  
  ###include ordinal thresholds
  if(!is.null(ctmodelobj$ordNames[1])){
    ordNames<-ctmodelobj$ordNames
    ordLevels<-ctmodelobj$ordLevels
    
    nOrdinal <- length(ordNames)
    ordvars <- ordNames
    
    if(objective!='Kalman' && objective != 'Kalmanmx') ordvars <- paste0(rep(ordNames,each=Tpoints),'_T',rep(0:(Tpoints-1),nOrdinal))
    if(objective=='Kalman' || objective == 'Kalmanmx') Tpoints <- 1
    maxthresholds <- max(unlist(lapply(ordLevels,length)))-1
    nthresholds <- (unlist(lapply(ordLevels,length)))-1
    tb <- matrix(NA,maxthresholds,length(ordvars))
    colnames(tb) <- ordvars
    tfree <- tb
    tfree[,] <- FALSE
    tvalues <- tb
    tlabels <- tb
    tlbound <- tb
    tubound <- tb
    
    dat <- as.data.frame(model@data$observed)
    
    
    for(vari in 1:length(ordNames)){
      tfree[1:nthresholds[vari],(1+(Tpoints*(vari-1))):(Tpoints*vari)] <- TRUE
      tvalues[1:nthresholds[vari],(1+(Tpoints*(vari-1))):(Tpoints*vari)] <- ordLevels[[vari]][-(1+nthresholds[vari])]+.5
      tlabels[1:nthresholds[vari],(1+(Tpoints*(vari-1))):(Tpoints*vari)] <- paste0('thresh_',ordNames[vari],'_',1:nthresholds[vari])
      tlbound[1:nthresholds[vari],(1+(Tpoints*(vari-1))):(Tpoints*vari)] <- ordLevels[[vari]][-(1+nthresholds[vari])]
      tubound[1:nthresholds[vari],(1+(Tpoints*(vari-1))):(Tpoints*vari)] <- ordLevels[[vari]][-1]
      dat[,ordvars[(1+(Tpoints*(vari-1))):(Tpoints*vari)]] <- mxFactor(dat[,ordvars[(1+(Tpoints*(vari-1))):(Tpoints*vari)]],levels = ordLevels[[vari]]) 
    }
    # browser()
    thresh <- mxMatrix(name='thresh',
      ncol=length(ordvars),
      nrow=maxthresholds,
      dimnames = list(paste0('row',1:maxthresholds),ordvars),
      free = tfree,
      values = tvalues,
      labels = tlabels,
      lbound = tlbound,
      ubound= tubound)
    
    model<-mxModel(
  model, 
  thresh,
  mxData(observed = dat, type = "raw"))
    
  }
  
  
  
  if(plotOptimization==TRUE){
    model<- OpenMx::mxOption(model, 'Always Checkpoint', 'Yes')
    model<- OpenMx::mxOption(model, 'Checkpoint Units', 'iterations')
    model<- OpenMx::mxOption(model, 'Checkpoint Count', 1)    
  }
  
  model<-mxOption(model,'RAM Inverse Optimization', 'No')
  model<-mxOption(model, 'Nudge zero starts', FALSE)
  # model<-mxOption(model, 'Gradient step size', 1e-4)
  
  ###fit model
  if(!is.null(omxStartValues)) model<-omxSetParameters(model,
    labels=names(omxStartValues)[names(omxStartValues) %in% names(omxGetParameters(model))],
    values=omxStartValues[names(omxStartValues) %in% names(omxGetParameters(model))],strict=FALSE)
  
  if(carefulFit==TRUE) {
    carefulFit<-FALSE
    
    mxobj<-carefulFit(model,traitExtension=traitExtension,
      manifestTraitvarExtension=manifestTraitvarExtension,weighting=carefulFitWeight) #fit with the penalised likelihood
    #         mxobj<-OpenMx::mxRun(model) #fit with the penalised likelihood
    
    # message(paste0('carefulFit penalisation:  ', mxEval(ctsem.penalties, mxobj,compute=TRUE),'\n'))
    newstarts <- try(OpenMx::omxGetParameters(mxobj)) #get the params
    if(showInits==TRUE) {
      message('Generated start values from carefulFit=TRUE')
      message(paste(names(newstarts), ": ", newstarts, "\n"))
    }
    
    if(class(newstarts)!="try-error" & !is.null(newstarts)) model<-OpenMx::omxSetParameters(model, 
      labels=names(newstarts),  values=newstarts,strict=FALSE) #set the params of it
    #     objective<-targetObjective #revert our objective to whatever was specified
    #     setobjective() #and set it
  }
  
  if(fit == FALSE) mxobj <- model #if we're not fitting the model, just return the unfitted openmx model
  
  if(fit == TRUE){ #but otherwise...  
    
    if(useOptimizer==TRUE) mxobj <- mxTryHard(model, initialTolerance=1e-14,
      # finetuneGradient=FALSE,
      initialGradientIterations=1,
      #initialGradientStepSize = 1e-6,
      showInits=showInits, checkHess=TRUE, greenOK=FALSE, 
      iterationSummary=iterationSummary, bestInitsOutput=FALSE, verbose=verbose,
      extraTries=retryattempts, loc=.5, scale=0.2, paste=FALSE)
    
    if(useOptimizer==FALSE) mxobj <- OpenMx::mxRun(model,useOptimizer=useOptimizer)
  }
  
  
  if(plotOptimization==TRUE){
    
    checkpoints<-utils::read.table(file='ctsem.omx', header=TRUE, sep='\t')
    if(carefulFit==TRUE) checkpoints<-rbind(checkpoints, NA, NA, NA, NA, NA, utils::read.table(file='ctsemCarefulfit.omx', header=TRUE, sep='\t'))
    mfrow<-graphics::par()$mfrow
    graphics::par(mfrow=c(3, 3))
    for(i in 6:ncol(checkpoints)) {
      graphics::plot(checkpoints[, i], main=colnames(checkpoints)[i])
    }
    graphics::par(mfrow=mfrow)
    deleteCheckpoints <- readline('Remove created checkpoint file, ctsem.omx? y/n \n')
    if(deleteCheckpoints=='y') file.remove(file='ctsem.omx')
  }
  
  OpenMx::mxOption(NULL, "Default optimizer", originaloptimizer) #reset optimizer
  
  out <- list(mxobj, ctmodelobj, ctfitargs, OpenMx::omxGetParameters(model), startValues) #roll unfitted and fitted model and ctmodelobj into one list item
  names(out) <- c("mxobj", "ctmodelobj", "ctfitargs", 'omxStartValues', 'ctStartValues')
  class(out) <- "ctsemFit" #and give it the class of a ctsemFit object
  
  return(out)
}
