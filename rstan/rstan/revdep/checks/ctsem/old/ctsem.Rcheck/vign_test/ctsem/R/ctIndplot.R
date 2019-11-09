#' ctIndplot
#' 
#' Convenience function to simply plot individuals trajectories from ctsem wide format data
#' @param datawide ctsem wide format data
#' @param n.subjects Number of subjects to randomly select for plotting, or character vector 'all'.
#' @param Tpoints Number of discrete time points per case in data structure
#' @param n.manifest Number of manifest variables in data structure
#' @param colourby set  plot colours by "subject" or "variable"
#' @param vars either 'all' or a numeric vector specifying which manifest variables to plot.
#' @param opacity Opacity of plot lines
#' @param varnames vector of variable names for legend (defaults to NULL)
#' @param xlab X axis label.
#' @param ylab Y axis label.
#' @param type character specifying plot type, as per usual base R plot commands. 
#' Defaults to 'b', both points and lines.
#' @param start Measurement occasion to start plotting from - defaults to T0.
#' @param legend Logical. Plot a legend?
#' @param legendposition Where to position the legend.
#' @param ... additional plotting parameters.
#' @param new logical. If TRUE, creates a new plot, otherwise overlays on current plot.
#' @param jittersd positive numeric indicating standard deviation of noise to add to observed
#' data for plotting purposes.
#' @examples 
#' 
#' data(ctExample1)
#' ctIndplot(ctExample1,n.subjects=1, n.manifest=2,Tpoints=6, colourby='variable')
#' 
#' @export
ctIndplot<-function(datawide,n.manifest,Tpoints,n.subjects='all',colourby="variable",
  vars='all',opacity=1,varnames=NULL,xlab='Time',ylab='Value',type='b',start=0,legend=TRUE,
  legendposition='topright',new=TRUE,jittersd=.05,...){

  if(n.subjects=='all') n.subjects=nrow(datawide)
  
  if(length(vars)==1 && vars=='all') vars<-1:n.manifest
  
  if(colourby=="variable") colourvector <- grDevices::rainbow(length(vars),alpha=opacity)
  if(colourby=="subject") colourvector <- grDevices::rainbow(n.subjects,alpha=opacity)
 
  
  ymin<-min(datawide[1:nrow(datawide),cseq(vars,n.manifest*Tpoints,n.manifest)],na.rm=T)
  ymax<-max(datawide[1:nrow(datawide),cseq(vars,n.manifest*Tpoints,n.manifest)],na.rm=T)
  
  # browser()
  individuals<-sample(1:nrow(datawide),n.subjects)
  times<-matrix(unlist(lapply(1:(Tpoints-1),function(x){
    apply(datawide[individuals,,drop=FALSE][,paste0('dT',1:x),drop=FALSE],1,sum,na.rm=T)
  })),ncol=(Tpoints-1))
  
  if(new==TRUE) graphics::plot(NA,ylim=c(ymin,ymax),xlim=c(start,max(times)),
    ylab=ylab,xlab=xlab,...)
  
 
  message(c('Plotting rows ',paste0(individuals,", ")))
  for(i in 1:n.subjects){

    for(j in 1:length(vars)){
      graphics::points(c(0,times[i,]),
        datawide[individuals[i],seq(vars[j],n.manifest*Tpoints,n.manifest)] +
          rnorm(length(datawide[individuals[i],seq(vars[j],n.manifest*Tpoints,n.manifest)]),0,jittersd),
        col=ifelse(colourby=="variable",colourvector[j],colourvector[i]),type=type,pch=j,lty=1,...) 
    }}
  
  if(is.null(varnames)) varnames <- paste0("var",vars) #set varnames for legend
  
  if(legend==TRUE){
  if(colourby=="variable") {
    graphics::legend(legendposition,varnames,pch=vars,col=colourvector,text.col=colourvector,bty="n")
  }
  if(colourby=="subject") {
    graphics::legend(legendposition,varnames,pch=vars,bty="n")
  }
  }
  
}
