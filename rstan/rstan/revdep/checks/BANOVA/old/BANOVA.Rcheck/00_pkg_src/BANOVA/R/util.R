assign_contrast <- function(data, # a data frame
                            contrast # must be a list of matrix
                            ){
  if (!is.null(contrast)){
    c_names = attr(contrast, 'names')
    for (c_name in c_names){
      # check the dim of each contrast
      if (is.null(dim(contrast[[c_name]]))){
        d = 1
      }else{
        d = dim(contrast[[c_name]])[2]
      }
      if (c_name %in% colnames(data))
        contrasts(data[[c_name]], d) = contrast[[c_name]]
    }
  }
  return(data)
}


# for table of predictions
assign_contrast_factor <- function(factor, # factor value
                                   name, # factor name
                                   contrast # must be a list of matrix
                                   ){
  if (!is.null(contrast)){
    c_names = attr(contrast, 'names')
    if (name %in% c_names){
      # check the dim of each contrast
      if (is.null(dim(contrast[[name]]))){
        d = 1
      }else{
        d = dim(contrast[[name]])[2]
      }
      contrasts(factor, d) = contrast[[name]]
    }
  }
  return(factor)
}

# for mediation/floodlight effects using effect coding
#
#

# create idmap new_id -> old_id
idmap <- function(old_id, new_id){
  # new_id is from 1,2,3,...
  # old_id is any character string
  res <- rep("", length(unique(new_id)))
  for (i in 1:length(new_id)){
    if (res[new_id[i]] == "")
      res[new_id[i]] = old_id[i]
  }
  return(res)
}