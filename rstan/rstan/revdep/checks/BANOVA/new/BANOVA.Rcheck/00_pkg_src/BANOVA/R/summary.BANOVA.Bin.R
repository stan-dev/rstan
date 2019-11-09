summary.BANOVA.Binomial <-
function(object, ...){
  sol <- list(anova.table = object$anova.table,
              coef.table = object$coef.tables$full_table,
              pvalue.table = object$pvalue.table, conv = object$conv, full_object = object, call = object$call)
  class(sol) <- "summary.BANOVA.Binomial"
  sol
}
