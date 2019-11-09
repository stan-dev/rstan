#' ctCI
#' Computes confidence intervals on specified parameters / matrices for already fitted ctsem fit object.
#' @param ctfitobj Already fit ctsem fit object (class: ctsemFit) to estimate confidence intervals for.
#' @param confidenceintervals character vector of matrices and or parameters for which to estimate 95\% confidence intervals for.
#' @param optimizer character vector. Defaults to NPSOL (recommended), but other optimizers available within OpenMx (e.g. 'CSOLNP') may be specified.
#' @param verbose Integer between 0 and 3 reflecting amount of output while calculating.
#' @return ctfitobj, with confidence intervals included.
#' @details Confidence intervals typically estimate more reliably using the proprietary NPSOL optimizer available within OpenMx only when
#' installing directly from OpenMx website. Use command " source('http://openmx.psyc.virginia.edu/getOpenMx.R') " to install OpenMx with NPSOL.
#' If estimating for a multigroup model, specify confidence intervals as normal, e.g. \code{confidenceintervals = c('DRIFT', 'diffusion_Y1_Y1')} .
#' The necessary group prefixes are added internally.
#' 
#' @examples
#' ## Examples set to 'dontrun' because they take longer than 5s.
#' \dontrun{
#' data("ctExample3")
#' model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
#'  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1), 
#'  MANIFESTMEANS = matrix(c(0, "manifestmean2", "manifestmean3"), nrow = 3, 
#'    ncol = 1))
#' fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = "Kalman",
#'  stationary = c("T0VAR"))
#' 
#' fit <- ctCI(fit, confidenceintervals = 'DRIFT')
#' 
#' summary(fit)$omxsummary$CI
#' }
#' @export
ctCI<-function(ctfitobj, confidenceintervals, optimizer='NPSOL', verbose=0){
  
  if(class(ctfitobj)!= 'ctsemFit' & class(ctfitobj)!= 'ctsemMultigroupFit') stop('Not a ctsem fit object!')
  
  if(class(ctfitobj)=='ctsemMultigroupFit') confidenceintervals<-unlist(lapply(confidenceintervals,function(x){
    paste0(rep(names(ctfitobj$mxobj$submodels),each=length(confidenceintervals)),'.', confidenceintervals)
  }))
  
  originalOptimizer<- mxOption(NULL, "Default optimizer")
  
  if(optimizer=='NPSOL' & imxHasNPSOL()==FALSE) {
    optimizer='CSOLNP'
    message("NPSOL optimizer not available - switching to CSOLNP - consider installing OpenMx using command:  source('http://openmx.psyc.virginia.edu/getOpenMx.R') ")
  }
    mxOption(NULL,'Default optimizer', optimizer)

  ctfitobj$mxobj <- OpenMx::mxModel(ctfitobj$mxobj, 
    mxCI(confidenceintervals, interval = 0.95, type = "both"),
  mxComputeSequence(list(
    mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
      maxMajorIter=3000),
      constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'SLSQP','ineq','none')),
    mxComputeNumericDeriv(), mxComputeStandardError(), 
    mxComputeReportDeriv()))
  )
  

  ctfitobj$mxobj<-mxRun(ctfitobj$mxobj,intervals=TRUE)
  ctfitobj$mxobj@compute<-NULL
  mxOption(NULL,'Default optimizer', originalOptimizer)
  
  return(ctfitobj)
}
