trace.plot <-
function(x, save = FALSE){
  if (class(x) %in% c('BANOVA', 'BANOVA.Binomial','BANOVA.Bernoulli',
                      'BANOVA.Normal', 'BANOVA.T', 'BANOVA.Multinomial','BANOVA.ordMutinomial')){
    if (!x$single_level){
      if (x$model_name == 'BANOVA.Multinomial')
        traceplot(x$samples_l2_param, colnames(x$dMatrice$X_full[[1]]), colnames(x$dMatrice$Z), save)
      else
        traceplot(x$samples_l2_param, colnames(x$dMatrice$X), colnames(x$dMatrice$Z), save)
    }else{
      if (x$model_name == 'BANOVA.Multinomial')
        traceplot(x$samples_l1_param, Z_names = colnames(x$dMatrice$X_full[[1]]), save = save)
      else
        traceplot(x$samples_l1_param, Z_names = colnames(x$dMatrice$X), save = save)
    }
  }
}
