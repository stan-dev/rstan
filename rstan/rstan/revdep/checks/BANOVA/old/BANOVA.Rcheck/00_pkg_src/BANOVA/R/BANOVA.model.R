#' extract BANOVA models
#' @param model_name String specifying BANOVA models (e.g. Normal)
#'
#' @return A BANOVA model (stan code).
#'
#' @examples
#' \dontrun{
#' 
#' }
#'
BANOVA.model <- function (model_name, 
                          single_level = F){
  if(model_name %in% c('Normal', 'Poisson', 'T', 'Bernoulli', 
                       'Binomial', 'ordMultinomial', 'Multinomial')){
    if(single_level){
        name <- paste("stan/single_",model_name, ".stan", sep = "")
    }else{
        name <- paste("stan/", model_name, "_Normal.stan", sep = "")
    }
    file_src <- system.file(name, package = 'BANOVA', mustWork = TRUE)
    #file_src <- paste("Windows/BANOVA_v1.1/BANOVA_R/inst/",name,sep = "")
    model_code = readChar(file_src, nchars=1000000)

  }else{
    stop(paste(model_name, " model is not supported currently!"))
  }
  sol <- list(model_code = model_code, model_name = model_name, single_level = single_level)
  class(sol) <- 'BANOVA.model'
  return(sol)
}

