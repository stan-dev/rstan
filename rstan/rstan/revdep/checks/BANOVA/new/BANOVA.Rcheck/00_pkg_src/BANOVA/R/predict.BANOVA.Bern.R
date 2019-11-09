predict.BANOVA.Bernoulli <-
function(object, newdata = NULL, ...){
  sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                       mf1 = object$mf1, mf2 = object$mf2, newdata, model = 'BernNormal', n_trials = object$num_trials)
  return(sol)
}
