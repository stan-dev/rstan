#' build BANOVA models
#' @param model_code BANOVA models in stan format
#'
#' @return A BANOVA stan model.
#'
#' @examples
#' \dontrun{
#' 
#' }
#'
BANOVA.build <- function(BANOVA_model){
  if(!is(BANOVA_model, 'BANOVA.model')) stop('BANOVA_model must be a BANOVA.model object, use the BANOVA.model function to create a model first!')
  cat('Compiling...\n')
  stan_c <- stanc(model_code = BANOVA_model$model_code, model_name = BANOVA_model$model_name)
  utils::capture.output(stanmodel <- rstan::stan_model(stanc_ret = stan_c,save_dso = T), type = "output")
  cat('Compiled successfully\n')
  sol <- list(stanmodel = stanmodel, model_name = BANOVA_model$model_name, single_level = BANOVA_model$single_level)
  class(sol) <- 'BANOVA.build'
  return(sol)
}