print.ancova.effect <-
  function(x, ...){
    cat("\nTable of sum of squares:\n")
    print(x$ancova_table)
    cat("\n")
    cat("Table of effect sizes (95% credible interval):\n")
    print(x$effect_table)
  }