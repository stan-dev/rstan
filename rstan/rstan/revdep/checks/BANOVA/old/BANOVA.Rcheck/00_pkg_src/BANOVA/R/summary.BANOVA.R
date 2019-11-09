summary.BANOVA <- function(object, ...){
  if(object$model_name == "BANOVA.ordMultinomial"){
    res <- list(anova.table = object$anova.table,
                coef.table = rbind(object$coef.tables$full_table, object$coef.tables$cutp_table),
                pvalue.table = object$pvalue.table, conv = object$conv, R2 = object$R2, full_object = object, 
                single_level = object$single_level, model_name = object$model_name, call = object$call)
  }else{
    res <- list(anova.table = object$anova.table,
                coef.table = object$coef.tables$full_table,
                pvalue.table = object$pvalue.table, conv = object$conv, R2 = object$R2, full_object = object, 
                single_level = object$single_level, model_name = object$model_name, call = object$call)
  }
  class(res) <- "summary.BANOVA"
  res
}