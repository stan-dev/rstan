pkgname <- "banocc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('banocc')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("banocc_model")
### * banocc_model

flush(stderr()); flush(stdout())

### Name: banocc_model
### Title: The stan model used in the Bayesian fit
### Aliases: banocc_model
### Keywords: datasets

### ** Examples

data(compositions_null)
## Not run: 
##D   compiled_banocc_model <- rstan::stan_model(model_code = banocc_model)
## End(Not run)




cleanEx()
nameEx("get_banocc_output")
### * get_banocc_output

flush(stderr()); flush(stdout())

### Name: get_banocc_output
### Title: Takes a model fit from BAnOCC, evaluates convergence and
###   generates appropriate convergence metrics and inference
### Aliases: get_banocc_output

### ** Examples

data(compositions_null)
  ## Not run: 
##D     compiled_banocc_model <- rstan::stan_model(model_code=banocc_model)
##D     b_fit <- run_banocc(C=compositions_null,
##D                             compiled_banocc_model=compiled_banocc_model)
##D     b_output <- get_banocc_output(banoccfit=b_fit)
##D   
## End(Not run)




cleanEx()
nameEx("run_banocc")
### * run_banocc

flush(stderr()); flush(stdout())

### Name: run_banocc
### Title: Runs BAnOCC to fit the model and generate appropriate
###   convergence metrics and inference.
### Aliases: run_banocc

### ** Examples

  data(compositions_null)
  ## Not run: 
##D     compiled_banocc_model <- rstan::stan_model(model_code=banocc_model)
##D     b_stanfit <- run_banocc(C=compositions_null,
##D                             compiled_banocc_model=compiled_banocc_model)
##D   
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
