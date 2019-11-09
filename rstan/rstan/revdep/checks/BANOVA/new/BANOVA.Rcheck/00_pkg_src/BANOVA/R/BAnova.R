BAnova <-
function(x){
  if (length(x$anova.table) >2){
    for (i in 1:length(x$anova.table)){
      cat('\nChoice: ', i, '\n')
      print(x$anova.table[[i]])
    }
  }else{
    print(x$anova.table)
  }
}
