predict.BANOVA.ordMultinomial <-
function(object, newdata = NULL,...){
  if(object$single_level){
    sol <- predict.means(object$samples_l1_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                         mf1 = object$mf1, mf2 = object$mf2, newdata, samples_cutp_param = object$samples_cutp_param, model = 'MultinomialordNormal')
  }else{
    sol <- predict.means(object$samples_l2_param, object$data, object$dMatrice$X, object$dMatrice$Z_full, 
                         mf1 = object$mf1, mf2 = object$mf2, newdata, samples_cutp_param = object$samples_cutp_param, model = 'MultinomialordNormal')
  }
  return(sol)
}
