predict.BANOVA <- function(object, newdata = NULL, Xsamples = NULL, Zsamples = NULL,...){
  if(object$single_level){
    if (object$model_name %in% c("BANOVA.Normal", "BANOVA.T")){
      sol <- predict.means(object$samples_l1_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'NormalNormal')
    }else if(object$model_name == "BANOVA.Poisson"){
      sol <- predict.means(object$samples_l1_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'PoissonNormal')
    }else if(object$model_name %in% c("BANOVA.Binomial", "BANOVA.Bernoulli")){
      sol <- predict.means(object$samples_l1_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'BernNormal', n_trials = object$num_trials)
    }else if(object$model_name == "BANOVA.ordMultinomial"){
      sol <- predict.means(object$samples_l1_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, samples_cutp_param = object$samples_cutp_param, model = 'MultinomialordNormal')
    }else if(object$model_name == "BANOVA.Multinomial"){
      sol <- multi.predict.means(object$samples_l1_param, object$dataX, object$dataZ, object$dMatrice$X_full, object$dMatrice$X_original_choice, object$dMatrice$Z_full, 
                                 mf1 = object$mf1, mf2 = object$mf2, Xsamples, Zsamples)
    }
  }else{
    if (object$model_name %in% c("BANOVA.Normal", "BANOVA.T")){
      sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'NormalNormal')
    }else if(object$model_name == "BANOVA.Poisson"){
      sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'PoissonNormal', l2_sd = object$samples_l2_sigma_param)
    }else if(object$model_name %in% c("BANOVA.Binomial", "BANOVA.Bernoulli")){
      sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'BernNormal', n_trials = object$num_trials)
    }else if(object$model_name == "BANOVA.ordMultinomial"){
      sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                           mf1 = object$mf1, mf2 = object$mf2, newdata, samples_cutp_param = object$samples_cutp_param, model = 'MultinomialordNormal')
    }else if(object$model_name == "BANOVA.Multinomial"){
      sol <- multi.predict.means(object$samples_l2_param, object$dataX, object$dataZ, object$dMatrice$X_full, object$dMatrice$X_original_choice, object$dMatrice$Z_full, 
                                 mf1 = object$mf1, mf2 = object$mf2, Xsamples, Zsamples)
    }
  }
  return(sol)
}