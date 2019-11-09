pkgname <- "CopulaDTA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CopulaDTA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cdtamodel")
### * cdtamodel

flush(stderr()); flush(stdout())

### Name: cdtamodel
### Title: Specify the copula based bivariate beta-binomial distribution to
###   fit to the diagnostic data.
### Aliases: cdtamodel

### ** Examples

data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))



cleanEx()
nameEx("fit.cdtamodel")
### * fit.cdtamodel

flush(stderr()); flush(stdout())

### Name: fit.cdtamodel
### Title: Fit copula based bivariate beta-binomial distribution to
###   diagnostic data.
### Aliases: fit.cdtamodel

### ** Examples

data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
## Not run: 
##D fit1 <- fit(model1,
##D                 SID='ID',
##D                 data=telomerase,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
##D 
##D 
##D fit2 <- fit(model2,
##D                 SID='StudyID',
##D                 data=ascus,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
## End(Not run)




cleanEx()
nameEx("forestplot.cdtafit")
### * forestplot.cdtafit

flush(stderr()); flush(stdout())

### Name: forestplot.cdtafit
### Title: Produce forest plots for categorical covariates.
### Aliases: forestplot.cdtafit

### ** Examples

data(telomerase)
model1 =  cdtamodel(copula = 'fgm')

model2 = cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 =  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
## Not run: 
##D fit1 <- fit(model1,
##D                 SID='ID',
##D                 data=telomerase,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
##D 
##D plot(fit1)
## End(Not run)



cleanEx()
nameEx("print.cdtafit")
### * print.cdtafit

flush(stderr()); flush(stdout())

### Name: print.cdtafit
### Title: Print a summary of the fitted model.
### Aliases: print.cdtafit

### ** Examples

data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
## Not run: 
##D 
##D fit1 <- fit(model1,
##D                 SID='ID',
##D                 data=telomerase,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
##D 
##D print(fit1)
##D 
## End(Not run)



cleanEx()
nameEx("summary.cdtafit")
### * summary.cdtafit

flush(stderr()); flush(stdout())

### Name: summary.cdtafit
### Title: Function to generate a summary a cdtafit object.
### Aliases: summary.cdtafit

### ** Examples

data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
## Not run: 
##D 
##D fit1 <- fit(model1,
##D                 SID='ID',
##D                 data=telomerase,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
##D 
##D ss <- summary(fit1)
##D 
## End(Not run)



cleanEx()
nameEx("traceplot.cdtafit")
### * traceplot.cdtafit

flush(stderr()); flush(stdout())

### Name: traceplot.cdtafit
### Title: Trace plot using ggplot2.
### Aliases: traceplot.cdtafit

### ** Examples

data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
## Not run: 
##D fit1 <- fit(model1,
##D                 SID='ID',
##D                 data=telomerase,
##D                 iter=2000,
##D                 warmup=1000,
##D                 thin=1,
##D                 seed=3)
##D 
##D traceplot(fit1)
##D 
##D traceplot(fit1) +
##D theme(axis.text.x = element_text(size=10, colour='black'),
##D       axis.text.y = element_text(size=10, colour='black'),
##D       axis.title.x = element_text(size=10, colour='black'),
##D       strip.text = element_text(size = 10, colour='black'),
##D       axis.title.y= element_text(size=10, angle=0, colour='black'),
##D       strip.text.y = element_text(size = 10, colour='black'),
##D       strip.text.x = element_text(size = 10, colour='black'),
##D       plot.background = element_rect(fill = "white", colour='white'),
##D       panel.grid.major = element_blank(),
##D       panel.background = element_blank(),
##D       strip.background = element_blank(),
##D       axis.line.x = element_line(color = 'black'),
##D       axis.line.y = element_line(color = 'black'))
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
