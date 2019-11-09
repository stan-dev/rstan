table.predictions <-
function(x){
  if (x$single_level){
    if(x$model_name == 'BANOVA.Normal' || x$model_name == 'BANOVA.T'){
      l2_values <- attr(x$dMatrice$X, 'varValues')
      l2_values[[1]] <- NULL  # remove y var
      l2_interactions <- attr(x$dMatrice$X, 'interactions')
      if (length(l2_interactions) > 0)
        for (i in 1: length(l2_interactions))
          l2_interactions[[i]] <- l2_interactions[[i]] - 1
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l1_param, X_names = colnames(x$dMatrice$Z), 
                                      X_assign = attr(x$dMatrice$Z, 'assign'),
                                      Z_names = colnames(x$dMatrice$X), 
                                      Z_assign = attr(x$dMatrice$X, 'assign'), 
                                      Z_classes = attr(x$dMatrice$X, 'dataClasses'),
                                      l2_values = l2_values,
                                      l2_interactions = l2_interactions, 
                                      l2_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                                      numeric_index_in_Z = attr(x$dMatrice$X, 'numeric_index'), 
                                      model = 'NormalNormal', contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.Poisson'){
      l2_values <- attr(x$dMatrice$X, 'varValues')
      l2_values[[1]] <- NULL  # remove y var
      l2_interactions <- attr(x$dMatrice$X, 'interactions')
      if (length(l2_interactions) > 0)
        for (i in 1: length(l2_interactions))
          l2_interactions[[i]] <- l2_interactions[[i]] - 1
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l1_param, X_names = colnames(x$dMatrice$Z), 
                                      X_assign = attr(x$dMatrice$Z, 'assign'),
                                      Z_names = colnames(x$dMatrice$X), 
                                      Z_assign = attr(x$dMatrice$X, 'assign'), 
                                      Z_classes = attr(x$dMatrice$X, 'dataClasses'),
                                      l2_values = l2_values,
                                      l2_interactions = l2_interactions, 
                                      l2_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                                      numeric_index_in_Z = attr(x$dMatrice$X, 'numeric_index'), 
                                      model = 'PoissonNormal', contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.Bernoulli' || x$model_name == 'BANOVA.Binomial' ){
      l2_values <- attr(x$dMatrice$X, 'varValues')
      l2_values[[1]] <- NULL  # remove y var
      l2_interactions <- attr(x$dMatrice$X, 'interactions')
      if (length(l2_interactions) > 0)
        for (i in 1: length(l2_interactions))
          l2_interactions[[i]] <- l2_interactions[[i]] - 1
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l1_param, X_names = colnames(x$dMatrice$Z), 
                                      X_assign = attr(x$dMatrice$Z, 'assign'),
                                      Z_names = colnames(x$dMatrice$X), 
                                      Z_assign = attr(x$dMatrice$X, 'assign'), 
                                      Z_classes = attr(x$dMatrice$X, 'dataClasses'),
                                      l2_values = l2_values,
                                      l2_interactions = l2_interactions, 
                                      l2_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                                      numeric_index_in_Z = attr(x$dMatrice$X, 'numeric_index'), 
                                      model = 'BernNormal',
                                      n_trials = x$num_trials, contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.ordMultinomial'){
      l2_values <- attr(x$dMatrice$X, 'varValues')
      l2_values[[1]] <- NULL  # remove y var
      l2_interactions <- attr(x$dMatrice$X, 'interactions')
      if (length(l2_interactions) > 0)
        for (i in 1: length(l2_interactions))
          l2_interactions[[i]] <- l2_interactions[[i]] - 1
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l1_param, X_names = colnames(x$dMatrice$Z), 
                                      X_assign = attr(x$dMatrice$Z, 'assign'),
                                      Z_names = colnames(x$dMatrice$X), 
                                      Z_assign = attr(x$dMatrice$X, 'assign'), 
                                      Z_classes = attr(x$dMatrice$X, 'dataClasses'),
                                      l2_values = l2_values,
                                      l2_interactions = l2_interactions, 
                                      l2_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                                      numeric_index_in_Z = attr(x$dMatrice$X, 'numeric_index'), 
                                      samples_cutp_param = x$samples_cutp_param,
                                      model = 'MultinomialordNormal', contrast = x$contrast)
      
    }else if (x$model_name == 'BANOVA.Multinomial'){
      # TODO: X_full[[1]] values must all be 1
      # factors are the same across all alternatives, 
      # which means different levels of the same factor has the same probability estimation. cancelled in prob/sum 
      # 
      sol_tables <- NULL
      # sol_tables <- multi.print.table.means (x$coef.tables$coeff_table, 
      #                                        n_choice = x$n_categories, 
      #                                        x$samples_l1_param, 
      #                                        X_names = colnames(x$dMatrice$X_full[[1]]), 
      #                                        X_assign = attr(x$dMatrice$X_full[[1]], 'assign'), 
      #                                        X_classes = attr(x$dMatrice$X_full[[1]], 'dataClasses'), 
      #                                        Z_names = colnames(x$dMatrice$Z), 
      #                                        Z_assign = attr(x$dMatrice$Z, 'assign'), 
      #                                        l1_values = attr(x$dMatrice$X_full[[1]], 'varValues'), 
      #                                        l1_interactions = attr(x$dMatrice$X_full[[1]], 'interactions'), 
      #                                        l1_interactions_index = attr(x$dMatrice$X_full[[1]], 'interactions_index'), 
      #                                        numeric_index_in_X = attr(x$dMatrice$X_full[[1]], 'numeric_index'),
      #                                        single_level = T)
     
    }
  }else{
    if(x$model_name == 'BANOVA.Normal' || x$model_name == 'BANOVA.T'){
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l2_param, colnames(x$dMatrice$X), X_assign = attr(x$dMatrice$X, 'assign'), 
                        X_classes = attr(x$dMatrice$X, 'dataClasses'), colnames(x$dMatrice$Z), Z_assign = attr(x$dMatrice$Z, 'assign'), 
                        Z_classes = attr(x$dMatrice$Z, 'dataClasses'), l1_values = attr(x$dMatrice$X, 'varValues'), 
                        l1_interactions = attr(x$dMatrice$X, 'interactions'), l1_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                        l2_values = attr(x$dMatrice$Z, 'varValues'), l2_interactions = attr(x$dMatrice$Z, 'interactions'), 
                        l2_interactions_index = attr(x$dMatrice$Z, 'interactions_index'), numeric_index_in_X = attr(x$dMatrice$X, 'numeric_index'),
                        numeric_index_in_Z = attr(x$dMatrice$Z, 'numeric_index'), model = 'NormalNormal', contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.Poisson'){
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l2_param, colnames(x$dMatrice$X), X_assign = attr(x$dMatrice$X, 'assign'), 
                        X_classes = attr(x$dMatrice$X, 'dataClasses'), colnames(x$dMatrice$Z), Z_assign = attr(x$dMatrice$Z, 'assign'), 
                        Z_classes = attr(x$dMatrice$Z, 'dataClasses'), l1_values = attr(x$dMatrice$X, 'varValues'), 
                        l1_interactions = attr(x$dMatrice$X, 'interactions'), l1_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                        l2_values = attr(x$dMatrice$Z, 'varValues'), l2_interactions = attr(x$dMatrice$Z, 'interactions'), 
                        l2_interactions_index = attr(x$dMatrice$Z, 'interactions_index'), numeric_index_in_X = attr(x$dMatrice$X, 'numeric_index'),
                        numeric_index_in_Z = attr(x$dMatrice$Z, 'numeric_index'), model = 'PoissonNormal', l2_sd = x$samples_l2_sigma_param, contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.Bernoulli' || x$model_name == 'BANOVA.Binomial' ){
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l2_param, colnames(x$dMatrice$X), X_assign = attr(x$dMatrice$X, 'assign'), 
                        X_classes = attr(x$dMatrice$X, 'dataClasses'), colnames(x$dMatrice$Z), Z_assign = attr(x$dMatrice$Z, 'assign'), 
                        Z_classes = attr(x$dMatrice$Z, 'dataClasses'), l1_values = attr(x$dMatrice$X, 'varValues'), 
                        l1_interactions = attr(x$dMatrice$X, 'interactions'), l1_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                        l2_values = attr(x$dMatrice$Z, 'varValues'), l2_interactions = attr(x$dMatrice$Z, 'interactions'), 
                        l2_interactions_index = attr(x$dMatrice$Z, 'interactions_index'), numeric_index_in_X = attr(x$dMatrice$X, 'numeric_index'),
                        numeric_index_in_Z = attr(x$dMatrice$Z, 'numeric_index'), model = 'BernNormal', n_trials = x$num_trials, contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.ordMultinomial'){
      sol_tables <- print.table.means(x$coef.tables$coeff_table, x$samples_l2_param, colnames(x$dMatrice$X), X_assign = attr(x$dMatrice$X, 'assign'), 
                        X_classes = attr(x$dMatrice$X, 'dataClasses'), colnames(x$dMatrice$Z), Z_assign = attr(x$dMatrice$Z, 'assign'), 
                        Z_classes = attr(x$dMatrice$Z, 'dataClasses'), l1_values = attr(x$dMatrice$X, 'varValues'), 
                        l1_interactions = attr(x$dMatrice$X, 'interactions'), l1_interactions_index = attr(x$dMatrice$X, 'interactions_index'), 
                        l2_values = attr(x$dMatrice$Z, 'varValues'), l2_interactions = attr(x$dMatrice$Z, 'interactions'), 
                        l2_interactions_index = attr(x$dMatrice$Z, 'interactions_index'), numeric_index_in_X = attr(x$dMatrice$X, 'numeric_index'),
                        numeric_index_in_Z = attr(x$dMatrice$Z, 'numeric_index'), samples_cutp_param = x$samples_cutp_param, model = 'MultinomialordNormal', contrast = x$contrast)
    }else if (x$model_name == 'BANOVA.Multinomial'){
      sol_tables <- multi.print.table.means (x$coef.tables$coeff_table, n_choice = x$n_categories, x$samples_l2_param, colnames(x$dMatrice$X_full[[1]]), X_assign = attr(x$dMatrice$X_full[[1]], 'assign'), 
                               X_classes = attr(x$dMatrice$X_full[[1]], 'dataClasses'), colnames(x$dMatrice$Z), Z_assign = attr(x$dMatrice$Z, 'assign'), 
                               Z_classes = attr(x$dMatrice$Z, 'dataClasses'), l1_values = attr(x$dMatrice$X_full[[1]], 'varValues'), 
                               l1_interactions = attr(x$dMatrice$X_full[[1]], 'interactions'), l1_interactions_index = attr(x$dMatrice$X_full[[1]], 'interactions_index'), 
                               l2_values = attr(x$dMatrice$Z, 'varValues'), l2_interactions = attr(x$dMatrice$Z, 'interactions'), 
                               l2_interactions_index = attr(x$dMatrice$Z, 'interactions_index'), numeric_index_in_X = attr(x$dMatrice$X_full[[1]], 'numeric_index'),
                               numeric_index_in_Z = attr(x$dMatrice$Z, 'numeric_index'), contrast = x$contrast)
    }
  }
  invisible(sol_tables)
}
