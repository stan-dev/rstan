predict.BANOVA.Multinomial <-
function(object, Xsamples = NULL, Zsamples = NULL,...){
  sol <- multi.predict.means(object$samples_l2_param, object$dataX, object$dataZ, object$dMatrice$X_full, object$dMatrice$X_original_choice, object$dMatrice$Z_full, 
                       mf1 = object$mf1, mf2 = object$mf2, Xsamples, Zsamples)
  return(sol)
}
