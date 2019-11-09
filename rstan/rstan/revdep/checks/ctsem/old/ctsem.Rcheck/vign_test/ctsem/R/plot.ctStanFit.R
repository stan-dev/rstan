#' plot.ctStanFit
#'
#' Plots for ctStanFit objects
#'
#' @param x Fit object from \code{\link{ctStanFit}}.
#' @param types Vector of character strings defining which plots to create.
#' 'all' plots all possible types, including: 'regression', 'kalman', 
#' 'priorcheck', 'trace', 'density','intervals'. 
#' @param wait Logical. Pause between plots?
#' @param ... Arguments to pass through to the specific plot functions. Bewar of clashes
#' may occur if types='all'. For details see the specific functions generating each type of plot.
#' @details This function is just a wrapper calling the necessary functions for plotting - it 
#' may be simpler in many cases to access those directly. They are:
#'  \code{\link{ctStanDiscretePars}},\code{\link{ctKalman}},
#'  \code{\link{ctStanPlotPost}},\code{stan_trace},
#'  \code{stan_dens},\code{stan_plot}
#'  rstan offers many plotting possibilities not available here, to use that functionality
#'  one must simply call the relevant rstan plotting function. Use \code{x$stanfit} as the stan fit object
#'  (where x is the name of your ctStanFit object). Because a ctStanFit object has many 
#'  parameters, the additional argument \code{pars=ctStanParnames(x,'pop_')} is recommended.
#'  This denotes population means, but see \code{\link{ctStanParnames}} for
#'  other options.
#' @return Nothing. Generates plots.
#' @aliases ctStanPlot plot.ctStanFit
#' @method plot ctStanFit
#' @examples
#' \dontrun{
#' plot(ctstantestfit,types=c('regression','kalman','priorcheck'))
#' 
#' ### complete example
#' plot(ctstantestfit)
#' 
#' #### example plot using rstan functions
#' rstan::stan_trace(ctstantestfit$stanfit, 
#' pars=ctStanParnames(ctstantestfit,'pop_DRIFT'))
#' }
#' @export
plot.ctStanFit <- function(x, types='all',wait=TRUE,...){
  
  continue=TRUE
  
  waitf<-function(){
    out<-TRUE
    if(wait==TRUE && length(types) > 0){
      check <- readline("Input [s] to stop, or leave blank and press [return] for next plot.")
      if(check=='s' || check=='S') out<-FALSE
    }
    return(out)
  }
  
  if(types[1]=='all') types <- c('trace','regression', 'kalman', 'priorcheck','trace','density','intervals')
  
  if('regression' %in% types && continue){
    message('Plotting model implied regression coeffcients conditional on time interval using ctStanDiscretePars')
    
    ctStanDiscretePars(x, plot=TRUE,...)
    types=types[types!='regression']
    continue<-waitf()
  }
  
  if('kalman' %in% types && continue){
    message('Plotting expectations from Kalman filter using ctKalman')
    
    ctKalman(x, plot=TRUE,...)
    types=types[types!='kalman']
    continue<-waitf()
  }
  
  if('priorcheck' %in% types && continue){
    message('Plotting prior and posterior densities using ctStanPlotPost')
    
    ctStanPlotPost(x,...)
    types=types[types!='priorcheck']
    continue<-waitf()
  }
  
  
  if('trace' %in% types && continue){
    message('Plotting sampling traces using stan_trace')
    print(rstan::stan_trace(x$stanfit,ctStanParnames(x,'pop_'),...))
    continue<-waitf()
    
    if(continue) p<-try(rstan::stan_trace(x$stanfit,ctStanParnames(x,'popsd'),...),silent=TRUE)
     if(class(p)[1]!='try-error') {
       print(p)
       continue<-waitf()
     }
    
    if(continue)  {
      p<-try(rstan::stan_trace(x$stanfit,ctStanParnames(x,'tipred_'),...),silent=TRUE)
      types=types[types!='trace']
    }
    if(class(p)[1]!='try-error') {
      print(p)
      continue<-waitf()
    }
  }
  
  if('density' %in% types && continue){
    message('Plotting posterior density estimates using stan_dens')
    print(rstan::stan_dens(x$stanfit,ctStanParnames(x,'pop_'),...))
    continue<-waitf()
    
    if(continue)  p=try(rstan::stan_dens(x$stanfit,ctStanParnames(x,'popsd'),...),silent=TRUE)
    if(class(p)[1]!='try-error') {
      print(p)
      continue<-waitf()
    }
    
    # if(continue)  {
    #   p= try(rstan::stan_dens(x$stanfit,ctStanParnames(x,'tipred_'),...),silent=TRUE)
    #   types=types[types!='density']
    # }
    # if(class(p)[1]!='try-error') {
    #   print(p)
    #   continue<-waitf()
    # }
  }
  
  if('intervals' %in% types && continue){
    message('Plotting posterior intervals and point estimates using stan_plot')
    print(rstan::stan_plot(x$stanfit,ctStanParnames(x,'pop_'),...))
    continue<-waitf()
    
    if(continue)  p=try(rstan::stan_plot(x$stanfit,ctStanParnames(x,'popsd'),...),silent=TRUE)
    if(class(p)[1]!='try-error') {
      print(p)
      continue<-waitf()
    }
    # 
    # if(continue)  {
    #   p=try(rstan::stan_plot(x$stanfit,ctStanParnames(x,'tipred_'),...),silent=TRUE)
    #   types=types[types!='intervals']
    # }
    # if(class(p)[1]!='try-error') {
    #   print(p)
    #   continue<-waitf()
    # }
  }
}



