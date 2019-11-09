pkgname <- "bmlm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bmlm')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("isolate")
### * isolate

flush(stderr()); flush(stdout())

### Name: isolate
### Title: Create isolated within- (and optionally between-) person
###   variables.
### Aliases: isolate

### ** Examples

# Create within-person deviations of work stressors in BLch9.
data(BLch9)
BLch9 <- isolate(BLch9, by = "id", value = "fwkstrs")
head(BLch9)  # Now has new column for within-person work stressors.




cleanEx()
nameEx("mlm")
### * mlm

flush(stderr()); flush(stdout())

### Name: mlm
### Title: Estimate a multilevel mediation model
### Aliases: mlm

### ** Examples

## Not run: 
##D ## Run example from Bolger and Laurenceau (2013)
##D data(BLch9)
##D fit <- mlm(BLch9)
##D mlm_summary(fit)
##D 
##D ### With priors
##D Priors <- list(dy = 10, dm = 10, a = 2, b = 2, cp = 2,
##D                tau_dy = 1, tau_dm = 1, tau_a = 1, tau_b = 1, tau_cp = 1,
##D                lkj_shape = 2)
##D fit <- mlm(BLch9, priors = Priors)
## End(Not run)




cleanEx()
nameEx("mlm_path_plot")
### * mlm_path_plot

flush(stderr()); flush(stdout())

### Name: mlm_path_plot
### Title: Plot 'bmlm"s mediation model as a path diagram
### Aliases: mlm_path_plot

### ** Examples

# Draw a template path diagram of the mediation model
mlm_path_plot()




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
