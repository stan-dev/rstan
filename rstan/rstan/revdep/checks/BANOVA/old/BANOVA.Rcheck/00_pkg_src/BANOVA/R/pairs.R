# pairs functions
#' @param x a 'BANOVA' object
#'
#' @return plot pairs graph
#'
#' @examples
#' \dontrun{
#' 
#' }
pairs.BANOVA<- function(x, ...){
  if(x$model_name %in% c('BANOVA.Normal','BANOVA.T','BANOVA.Bernoulli','BANOVA.Binomial',
                     'BANOVA.Poisson','BANOVA.ordMultinomial','BANOVA.Multinomial')){
    pairs(x$stan_fit, ...)
  }else{
    stop(x$model_name, " is not supported by BANOVA.pairs!")
  }
}