print.BANOVA.mediation <- 
  function(x, ...){
    #cat('Direct effect of', x$xvar,':\n')
    #for (i in 1:length(x$dir_effects))
      #print(noquote(x$dir_effects[[i]]), row.names = F)
    if (x$individual){
      cat('Indirect effect of', x$xvar,':\n')
      for (i in 1:length(x$indir_effects)){
        print(noquote(x$individual_indirect[[i]]), row.names = F)
      }
    }else{
      cat('Indirect effect of', x$xvar,':\n')
      for (i in 1:length(x$indir_effects)){
        print(noquote(x$indir_effects[[i]]), row.names = F)
        cat('effect size: ', x$effect_size[[i]], '\n')
      }
    }
  }