packageStartupMessage('Compiling models which will take a while...')
lib_dir <- file.path(R_PACKAGE_DIR, paste('libs', R_ARCH, sep = ""))
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
packageStartupMessage(paste('Writing models to:', lib_dir))
#packageStartupMessage('No models compiled currently!')
BANOVA.model <- function (model_name, 
                          single_level = F){
  if(model_name %in% c('Normal', 'Poisson', 'T', 'Bernoulli', 
                       'Binomial', 'ordMultinomial', 'Multinomial')){
    if(single_level){
      name <- paste("single_",model_name, ".stan", sep = "")
    }else{
      name <- paste(model_name, "Normal.stan", sep = "_")
    }
    file_src <- file.path(R_PACKAGE_SOURCE, 'inst', 'stan', name)
    #file_src <- system.file(name, package = 'BANOVA', mustWork = TRUE)
    #file_src <- paste("BANOVA_v9/BANOVA_R/inst/stan/",name,sep = "")
    model_code = readChar(file_src, nchars=1000000)
    
  }else{
    stop(paste(model_name, " model is not supported currently!"))
  }
  sol <- list(model_code = model_code, model_name = model_name, single_level = single_level)
  class(sol) <- 'BANOVA.model'
  return(sol)
}

BANOVA.build <- function(BANOVA_model){
  if(!is(BANOVA_model, 'BANOVA.model')) stop('BANOVA_model must be a BANOVA.model object, use the BANOVA.model function to create a model first!')
  cat('Compiling the', BANOVA_model$model_name, 'model...\n')
  suppressMessages({
    stan_c <- rstan::stanc(model_code = BANOVA_model$model_code, model_name = BANOVA_model$model_name)
    utils::capture.output(stanmodel <- rstan::stan_model(stanc_ret = stan_c,save_dso = T), type = "output")
  })
  cat('Compiled successfully\n')
  sol <- list(stanmodel = stanmodel, model_name = BANOVA_model$model_name, single_level = BANOVA_model$single_level)
  class(sol) <- 'BANOVA.build'
  return(sol)
}

model_Normal_2 <- BANOVA.model('Normal')
Normal_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_Normal_2)

model_T_2 <- BANOVA.model('T', single_level = F)
T_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_T_2)

model_P_2 <- BANOVA.model('Poisson', single_level = F)
Poisson_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_P_2)

model_Bern_2 <- BANOVA.model('Bernoulli', single_level = F)
Bernoulli_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_Bern_2)

model_Bin_2 <- BANOVA.model('Binomial', single_level = F)
Binomial_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_Bin_2)

model_ordMulti_2 <- BANOVA.model('ordMultinomial', single_level = F)
ordMultinomial_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_ordMulti_2 )

model_Multi_2 <- BANOVA.model('Multinomial', single_level = F)
Multinomial_Normal_stanmodel <- BANOVA.build(BANOVA_model = model_Multi_2)

model_l2 <- file.path(lib_dir, 'BANOVA.RData')
save(Normal_Normal_stanmodel, Binomial_Normal_stanmodel, Bernoulli_Normal_stanmodel,
     T_Normal_stanmodel, Poisson_Normal_stanmodel, ordMultinomial_Normal_stanmodel,
     Multinomial_Normal_stanmodel, file = model_l2)

model_Normal_1 <- BANOVA.model('Normal', single_level = T)
Normal_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_Normal_1)

model_T_1 <- BANOVA.model('T', single_level = T)
T_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_T_1)

model_P_1 <- BANOVA.model('Poisson', single_level = T)
Poisson_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_P_1)

model_Bern_1 <- BANOVA.model('Bernoulli', single_level = T)
Bernoulli_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_Bern_1)

model_Bin_1 <- BANOVA.model('Binomial', single_level = T)
Binomial_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_Bin_1)

model_ordMulti_1 <- BANOVA.model('ordMultinomial', single_level = T)
ordMultinomial_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_ordMulti_1 )

model_Multi_1 <- BANOVA.model('Multinomial', single_level = T)
Multinomial_Normal_stanmodel_1 <- BANOVA.build(BANOVA_model = model_Multi_1)

model_l1 <- file.path(lib_dir, 'single_BANOVA.RData')
save(Normal_Normal_stanmodel_1, Binomial_Normal_stanmodel_1, Bernoulli_Normal_stanmodel_1,
     T_Normal_stanmodel_1, Poisson_Normal_stanmodel_1, ordMultinomial_Normal_stanmodel_1,
     Multinomial_Normal_stanmodel_1, file = model_l1)

packageStartupMessage('Models compiled successfully!')
