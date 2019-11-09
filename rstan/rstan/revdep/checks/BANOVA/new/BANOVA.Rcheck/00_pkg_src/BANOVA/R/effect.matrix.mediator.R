# find design matrix of coefficients of a mediator (including interactions which contains the mediator, using summation)
# the mediator must be included in each interaction
effect.matrix.mediator <- function (interaction_factors=list(), 
                                    mediator=NA, 
                                    matrix_formula=NULL, 
                                    xvar=NA, 
                                    xvar_include = FALSE, 
                                    intercept_include = FALSE, 
                                    flood_values = list(),
                                    contrast = NULL){
  # generate the effect matrix for each interaction, numerical covariates excluded, TODO: may have various levels for numeric variables
  # Args:
  #     interaction_factors       : a list of values of factors in the interaction to generate the effect matrix
  #     assign        : Not used for the new method, to remove
  #     index_factor  : Not used, to remove
  #     numeric_index : Not used, to remove
  # Returns:
  #     an effect matrix
  #
  tmp_contrasts <- getOption("contrasts")
  options(contrasts = rep("contr.sum",2))
  if (length(interaction_factors) > 0){
    num_inter = length(interaction_factors)
    levels_inter = list()
    names_inter = list()
    for (i in 1:num_inter){
      names_inter[[i]] = attr(interaction_factors[[i]], 'var_names')
      if (is.factor(interaction_factors[[i]])){ 
        levels_inter[[names_inter[[i]]]] = levels(interaction_factors[[i]])
      }else{
        # for numeric or integer var, use 1 instead, e.g. the mediator
        if(names_inter[[i]] %in% c(mediator, xvar)){
          levels_inter[[names_inter[[i]]]] = 1
        }else{
          if (names_inter[[i]] %in% names(flood_values)){
            if (length(names(flood_values)) != length(unique(names(flood_values)))) stop('Please provide unique flood light values!', call. = FALSE)
            if (length(flood_values[[names_inter[[i]]]]) > 1) stop('Please provide one value at a time for the other numeric variables!', call. = FALSE)
            levels_inter[[names_inter[[i]]]] = flood_values[[names_inter[[i]]]]
          }else{
            levels_inter[[names_inter[[i]]]] = 0 #mean(interaction_factors[[i]])
          }
        }
      }
    }
    factors_inter = expand.grid(levels_inter)
    level_factor <- factors_inter
    # remove the response variable
    if (!is.null(matrix_formula)){
      the_term <- terms(matrix_formula)
      if (attr(the_term, 'response') > 0){
        all_names <- as.character(attr(the_term, "variables"))[-1]
        response_name <- all_names[attr(the_term, 'response')]
        level_factor[response_name] <- NULL
      }
    }
    temp_fun <- function(x){
      if (length(unique(x)) > 1)
        as.factor(x)
      else
        as.numeric(as.character(x))
    }
    factors_inter_factor = as.data.frame(lapply(factors_inter, temp_fun))
    #### 1.1.2
    factors_inter_factor <- assign_contrast(factors_inter_factor, contrast)
    ####
    formula_inter = paste(names_inter, collapse = '*')
    if (is.null(matrix_formula))
      eval(parse(text = paste('model_frame <- model.frame(~',formula_inter,', data = factors_inter_factor)', sep='')))
    else
      model_frame <- model.frame(matrix_formula, data = factors_inter_factor)
      #eval(parse(text = paste('model_frame <- model.frame(',matrix_formula,', data = factors_inter_factor)', sep='')))

    effect_matrix <- model.matrix(formula(attr(model_frame, 'terms')), data = factors_inter_factor)
    #attr(effect_matrix, 'levels') = as.matrix(factors_inter)
    attr(effect_matrix, 'levels') = as.matrix(level_factor)
    if (!is.na(mediator)){
      # only calculate the summation of coefficients that related to the mediator exclude xvar
      if (!is.na(xvar) & (xvar %in% rownames(attr(attr(model_frame, 'terms'), 'factors')))){
        tmp_mtx <- attr(attr(model_frame, 'terms'), 'factors')
        tmp_ind_m <- which(rownames(tmp_mtx) == mediator)
        tmp_ind_x <- which(rownames(tmp_mtx) == xvar)
        if (!xvar_include)
          assign_selected <- which(tmp_mtx[tmp_ind_m, ] == 1 & tmp_mtx[tmp_ind_x, ] == 0) 
        else
          assign_selected <- which(tmp_mtx[tmp_ind_m, ] == 1 & tmp_mtx[tmp_ind_x, ] == 1) 
        #assign_selected <- which(attr(attr(model_frame, 'terms'), 'factors')[mediator, ] == 1 & attr(attr(model_frame, 'terms'), 'factors')[xvar, ] == 0) 
      }else{
        tmp_mtx <- attr(attr(model_frame, 'terms'), 'factors')
        tmp_ind <- which(rownames(tmp_mtx) == mediator)
        assign_selected <- which(tmp_mtx[tmp_ind, ] == 1)
        #assign_selected <- which(attr(attr(model_frame, 'terms'), 'factors')[mediator, ] == 1)
      }
      column_selected <- which(attr(effect_matrix, 'assign') %in% assign_selected)
      if (intercept_include) column_selected <- c(1,column_selected)
      effect_matrix_selected <- effect_matrix[ , column_selected, drop=FALSE]
      attr(effect_matrix_selected, 'levels') = as.matrix(level_factor)
    }else{
      if (!is.na(xvar) & (xvar %in% rownames(attr(attr(model_frame, 'terms'), 'factors')))){
        # exclude xvar
        tmp_mtx <- attr(attr(model_frame, 'terms'), 'factors')
        tmp_ind_x <- which(rownames(tmp_mtx) == xvar)
        assign_selected <- which(tmp_mtx[tmp_ind_x, ] == 0)
        column_selected <- which(attr(effect_matrix, 'assign') %in% assign_selected)
        if (length(column_selected) == 0){
          # intercept only
          effect_matrix_selected <- array(1, dim = c(1, 1), dimnames = list(NULL, '(Intercept)'))
          attr(effect_matrix_selected, 'levels') = as.matrix(array(1, dim = c(1, 1), dimnames = list(NULL, '(Intercept)')))
        }else{
          if (intercept_include) column_selected <- c(1,column_selected)
          effect_matrix_selected <- effect_matrix[ , column_selected, drop=FALSE]
          attr(effect_matrix_selected, 'levels') = as.matrix(level_factor)
        }
      }else{
        effect_matrix_selected <- effect_matrix
        attr(effect_matrix_selected, 'levels') = as.matrix(level_factor)
      }
    }
    
  }else{
    effect_matrix_selected = NA
  }
  options(contrasts = tmp_contrasts)
  return(effect_matrix_selected)
}
