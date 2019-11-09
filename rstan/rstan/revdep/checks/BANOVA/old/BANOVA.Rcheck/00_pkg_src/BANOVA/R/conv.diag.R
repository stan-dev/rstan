conv.diag <-
function(x){
  if (class(x) %in% c('BANOVA', 'BANOVA.Binomial','BANOVA.Bernoulli',
                      'BANOVA.Normal', 'BANOVA.T', 'BANOVA.Multinomial','BANOVA.ordMutinomial')){
    if (!x$single_level){
      if (x$model_name == 'BANOVA.Multinomial')
        sol <- conv.geweke.heidel(x$samples_l2_param, colnames(x$dMatrice$X_full[[1]]), colnames(x$dMatrice$Z))
      else
        sol <- conv.geweke.heidel(x$samples_l2_param, colnames(x$dMatrice$X), colnames(x$dMatrice$Z))
    }else{
      if (x$model_name == 'BANOVA.Multinomial')
        sol <- conv.geweke.heidel(x$samples_l1_param, colnames(x$dMatrice$Z), colnames(x$dMatrice$X_full[[1]]))
      else
        sol <- conv.geweke.heidel(x$samples_l1_param, colnames(x$dMatrice$Z), colnames(x$dMatrice$X))
    }
  }
  class(sol) <- 'conv.diag'
  sol
}
