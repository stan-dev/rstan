pkgname <- "BMSC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BMSC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("addInteractionToVars")
### * addInteractionToVars

flush(stderr()); flush(stdout())

### Name: addInteractionToVars
### Title: Add interactions of a specific order to a vector of variables
### Aliases: addInteractionToVars

### ** Examples

BMSC:::addInteractionToVars(3, c("x1", "x2", "x3"))



cleanEx()
nameEx("addPowToVars")
### * addPowToVars

flush(stderr()); flush(stdout())

### Name: addPowToVars
### Title: Add exponent to a vector of variables
### Aliases: addPowToVars

### ** Examples

BMSC:::addPowToVars(c("x1", "x2"), 2)



cleanEx()
nameEx("constrSelEst")
### * constrSelEst

flush(stderr()); flush(stdout())

### Name: constrSelEst
### Title: Model selection algorithm for constrained estimation
### Aliases: constrSelEst

### ** Examples

## Not run: 
##D suppressWarnings(RNGversion("3.5.0"))
##D set.seed(44)
##D n <- 80
##D x1 <- rnorm(n, sd = 1)
##D x2 <- rnorm(n, sd = 1)
##D x3 <- rnorm(n, sd = 1)
##D y <- 0.4 + 0.3 * x1 + 0.3 * x1 * x3 + 0.4 * x1 ^ 2 * x2 ^ 3 + rnorm(n, sd = 0.3)
##D yUncertainty <- rexp(n, 10) * 0.01
##D #optional (slow)
##D #xUncertainty <- data.frame(x3 = rep(0.1, n), x1 = rep(0.1, n), x2 = rep(1, n))
##D data <- data.frame(x1, x2, x3, y, yUncertainty)
##D models <- constrSelEst(y ~ x1 + x2 + x3, mustInclude = "x1", maxExponent = 3,
##D                        interactionDepth = 3, intercept = TRUE,
##D                        constraint_1 = TRUE, data = data,
##D                        yUncertainty = yUncertainty,
##D                        xUncertainty = NULL,
##D                        maxNumTerms = 10)
##D plotModelFit(models)
##D bestModel <- getBestModel(models, thresholdSE = 2)
##D print(bestModel)
## End(Not run)



cleanEx()
nameEx("createFormula")
### * createFormula

flush(stderr()); flush(stdout())

### Name: createFormula
### Title: Create a formula with interactions and polynomials up to a
###   desired order
### Aliases: createFormula

### ** Examples

createFormula("y ~ x1 + x2", 2, 3)
createFormula(as.formula("y ~ x1 + x2"), interactionDepth = 2)

carFormula <- createFormula("mpg ~ cyl + disp + drat", 2, 3)
summary(lm(carFormula, mtcars))



cleanEx()
nameEx("extractVarname")
### * extractVarname

flush(stderr()); flush(stdout())

### Name: extractVarname
### Title: Extract variable name from polynomial expression
### Aliases: extractVarname

### ** Examples

BMSC:::extractVarname(c("x1",
"I(x2^2)"))



cleanEx()
nameEx("makeInteractions")
### * makeInteractions

flush(stderr()); flush(stdout())

### Name: makeInteractions
### Title: Add all interactions up to a desired order
### Aliases: makeInteractions

### ** Examples

BMSC:::makeInteractions(vars = c("x1", "x2",
"I(x1^2)", "I(x2^2)"), interactionDepth = 3)



cleanEx()
nameEx("makePoly")
### * makePoly

flush(stderr()); flush(stdout())

### Name: makePoly
### Title: Create polynomial of degree 'maxExponent' from variable names
### Aliases: makePoly

### ** Examples

BMSC:::makePoly(vars = c("x1", "x2"), maxExponent = 3)



cleanEx()
nameEx("sortAndPaste")
### * sortAndPaste

flush(stderr()); flush(stdout())

### Name: sortAndPaste
### Title: Sort a vector and collapse elements together using ":"
### Aliases: sortAndPaste

### ** Examples

BMSC:::sortAndPaste(c("var1", "var2"))



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
