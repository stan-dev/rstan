pkgname <- "BayesSenMC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BayesSenMC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bd_meta")
### * bd_meta

flush(stderr()); flush(stdout())

### Name: bd_meta
### Title: Meta-analysis data on Bipolar Disorder diagnosis accuracy
### Aliases: bd_meta
### Keywords: dataset

### ** Examples

data(bd_meta)



cleanEx()
nameEx("correctedOR")
### * correctedOR

flush(stderr()); flush(stdout())

### Name: correctedOR
### Title: Model without misclassification
### Aliases: correctedOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}

# 3 MCMC chains with 10000 iterations each



cleanEx()
nameEx("crudeOR")
### * crudeOR

flush(stderr()); flush(stdout())

### Name: crudeOR
### Title: Model with constant nondifferential misclassification
### Aliases: crudeOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}\



cleanEx()
nameEx("diffOR")
### * diffOR

flush(stderr()); flush(stdout())

### Name: diffOR
### Title: Model with differential misclassification
### Aliases: diffOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}



cleanEx()
nameEx("fixedCorrOR")
### * fixedCorrOR

flush(stderr()); flush(stdout())

### Name: fixedCorrOR
### Title: Model with nondifferential, correlated misclassification
### Aliases: fixedCorrOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}



cleanEx()
nameEx("logitOR")
### * logitOR

flush(stderr()); flush(stdout())

### Name: logitOR
### Title: Model with nondifferential, logit normal-distributed
###   misclassification
### Aliases: logitOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}



cleanEx()
nameEx("nlmeNDiff")
### * nlmeNDiff

flush(stderr()); flush(stdout())

### Name: nlmeNDiff
### Title: Non-differential Generalized Linear Mixed Effects Model
### Aliases: nlmeNDiff

### ** Examples

data(bd_meta)

mod <- nlmeNDiff(bd_meta, lower = 0)



cleanEx()
nameEx("paramEst")
### * paramEst

flush(stderr()); flush(stdout())

### Name: paramEst
### Title: Parameter estimates of the GLMM model
### Aliases: paramEst

### ** Examples

data(bd_meta)

mod <- nlmeNDiff(bd_meta, lower = 0) # see nlme_nondiff() for detailed example.
pList <- paramEst(mod)



cleanEx()
nameEx("plotOR")
### * plotOR

flush(stderr()); flush(stdout())

### Name: plotOR
### Title: Plot Model
### Aliases: plotOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}

library(ggplot2)



cleanEx()
nameEx("randCorrOR")
### * randCorrOR

flush(stderr()); flush(stdout())

### Name: randCorrOR
### Title: Model with nondifferential, randomly correlated
###   misclassification
### Aliases: randCorrOR

### ** Examples

# Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
# Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}



cleanEx()
nameEx("smoke_meta")
### * smoke_meta

flush(stderr()); flush(stdout())

### Name: smoke_meta
### Title: Meta-analysis data on self-reported smoking diagnosis accuracy
### Aliases: smoke_meta
### Keywords: dataset

### ** Examples

data(smoke_meta)



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
