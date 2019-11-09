tformshapes <- function(){
paste0('  
  if(transform==0) out = param * meanscale * multiplier +inneroffset + offset; \n
  if(transform==1) out = log(1+(exp(param * meanscale+inneroffset))) * multiplier + offset ; \n
  if(transform==2) out = exp(param * meanscale+inneroffset) * multiplier + offset; \n
  if(transform==3) out = inv_logit(param*meanscale+inneroffset) * multiplier + offset; \n
  if(transform==4) out = ((param*meanscale+inneroffset)^3)*multiplier + offset; \n
')
  }

tform <- function(param, transform, multiplier, meanscale, offset, inneroffset, extratforms=''){
  if(!is.na(suppressWarnings(as.integer(transform)))) eval(parse(text=paste0(tformshapes(),extratforms)))
  if(is.na(suppressWarnings(as.integer(transform)))) out <- eval(parse(text=transform))  
  return(out)
}

