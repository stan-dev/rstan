pkgname <- "conStruct"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('conStruct')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("conStruct")
### * conStruct

flush(stderr()); flush(stdout())

### Name: conStruct
### Title: Run a conStruct analysis.
### Aliases: conStruct

### ** Examples

# load example dataset
data(conStruct.data)

# run example spatial analysis with K=1
	#	
# for this example, make.figs and save.files
#	are set to FALSE, but most users will want them 
#	set to TRUE
my.run <- conStruct(spatial = TRUE,
		 			K = 1,
		 			freqs = conStruct.data$allele.frequencies,
		 			geoDist = conStruct.data$geoDist,
		 			coords = conStruct.data$coords,
		 			prefix = "test",
		 			n.chains = 1,
		 			n.iter = 1e3,
		 			make.figs = FALSE,
		 			save.files = FALSE)




cleanEx()
nameEx("make.admix.pie.plot")
### * make.admix.pie.plot

flush(stderr()); flush(stdout())

### Name: make.admix.pie.plot
### Title: Make admixture pie plot
### Aliases: make.admix.pie.plot

### ** Examples

## Don't show: 
	admix.props <- matrix(c(0.086,0.000,0.500,0.505,0.099,0.052,0.024,0.007,0.800,0.000,0.216,0.744,0.917,0.199,0.469,0.000,0.783,0.298,0.329,0.446,0.000,0.000,0.637,0.903,0.000,0.000,0.000,0.012,0.021,0.000,0.000,0.089,0.000,0.554,0.002,0.000,0.000,0.095,0.020,0.001,0.001,0.011,0.000,0.200,0.000,0.060,0.053,0.082,0.036,0.013,0.000,0.062,0.169,0.137,0.029,0.001,0.000,0.178,0.079,0.000,0.999,1.000,0.988,0.979,0.975,1.000,0.744,0.984,0.435,0.998,0.914,1.000,0.405,0.475,0.900,0.947,0.965,0.993,0.000,1.000,0.725,0.203,0.000,0.765,0.518,1.000,0.154,0.533,0.534,0.525,0.999,1.000,0.185,0.018,1.000,0.001,0.000,0.000,0.000,0.025,0.000,0.167,0.016,0.012,0.000),ncol=3)
	coords <- matrix(c(-126.38,-125.23,-126.97,-128.54,-126.95,-121.71,-126.79,-123.38,-137.88,-125.82,-122.94,-130.73,-123.08,-122.84,-128.58,-124.82,-129.75,-122.25,-122.32,-129.10,-125.28,-123.98,-133.35,-131.74,-124.16,-146.35,-94.63,-149.02,-111.50,-126.67,-133.77,-118.63,-115.78,-113.42,-135.33,52.40,49.84,54.66,54.65,51.69,49.44,52.82,50.05,59.52,51.34,45.81,56.81,44.71,50.24,54.14,51.04,56.68,52.98,54.04,55.34,50.64,50.23,58.76,57.30,50.54,64.90,56.35,63.87,56.92,65.23,68.38,54.75,60.80,50.82,60.70),ncol=2)
## End(Don't show)	
# make admixture pie plot
make.admix.pie.plot(admix.proportions = admix.props,coords = coords)




cleanEx()
nameEx("make.structure.plot")
### * make.structure.plot

flush(stderr()); flush(stdout())

### Name: make.structure.plot
### Title: Make STRUCTURE output plot
### Aliases: make.structure.plot

### ** Examples

## Don't show: 
	admix.props <- matrix(c(0.086,0.000,0.500,0.505,0.099,0.052,0.024,0.007,0.800,0.000,0.216,0.744,0.917,0.199,0.469,0.000,0.783,0.298,0.329,0.446,0.000,0.000,0.637,0.903,0.000,0.000,0.000,0.012,0.021,0.000,0.000,0.089,0.000,0.554,0.002,0.000,0.000,0.095,0.020,0.001,0.001,0.011,0.000,0.200,0.000,0.060,0.053,0.082,0.036,0.013,0.000,0.062,0.169,0.137,0.029,0.001,0.000,0.178,0.079,0.000,0.999,1.000,0.988,0.979,0.975,1.000,0.744,0.984,0.435,0.998,0.914,1.000,0.405,0.475,0.900,0.947,0.965,0.993,0.000,1.000,0.725,0.203,0.000,0.765,0.518,1.000,0.154,0.533,0.534,0.525,0.999,1.000,0.185,0.018,1.000,0.001,0.000,0.000,0.000,0.025,0.000,0.167,0.016,0.012,0.000),ncol=3)
## End(Don't show)	
# make STRUCTURE-style plot
	make.structure.plot(admix.proportions = admix.props)

# make STRUCTURE-style plot, sorted by membership in layer 1
make.structure.plot(admix.proportions = admix.props,sort.by=1) 




cleanEx()
nameEx("match.layers.x.runs")
### * match.layers.x.runs

flush(stderr()); flush(stdout())

### Name: match.layers.x.runs
### Title: Match layers up across independent conStruct runs
### Aliases: match.layers.x.runs

### ** Examples

## Don't show: 
	admix.props1 <- matrix(c(0.09,0.00,0.50,0.51,0.10,0.05,0.02,0.01,0.80,0.00,0.22,0.74,0.92,0.20,0.47,0.00,0.78,0.30,0.33,0.45,0.00,0.00,0.64,0.90,0.00,0.00,0.00,0.01,0.02,0.00,0.00,0.09,0.00,0.55,0.00,0.00,0.00,0.09,0.02,0.00,0.00,0.01,0.00,0.20,0.00,0.06,0.05,0.08,0.04,0.01,0.00,0.06,0.17,0.14,0.03,0.00,0.00,0.18,0.08,0.00,1.00,1.00,0.99,0.98,0.98,1.00,0.74,0.98,0.43,1.00,0.91,1.00,0.41,0.47,0.90,0.95,0.96,0.99,0.00,1.00,0.72,0.20,0.00,0.77,0.52,1.00,0.15,0.53,0.53,0.53,1.00,1.00,0.18,0.02,1.00,0.00,0.00,0.00,0.00,0.02,0.00,0.17,0.02,0.01,0.00),ncol=3)
	admix.props2 <- matrix(c(0.36,0.35,0.42,0.38,0.35,0.35,0.36,0.35,0.48,0.36,0.39,0.39,0.40,0.36,0.36,0.35,0.40,0.46,0.45,0.38,0.34,0.35,0.47,0.40,0.35,1.00,1.00,0.99,0.99,0.98,1.00,0.84,0.99,0.63,1.00,0.32,0.35,0.24,0.24,0.33,0.34,0.33,0.35,0.15,0.32,0.32,0.10,0.30,0.33,0.27,0.36,0.13,0.26,0.27,0.22,0.36,0.35,0.14,0.11,0.35,0.00,0.00,0.00,0.01,0.01,0.00,0.07,0.00,0.18,0.00,0.32,0.30,0.34,0.38,0.31,0.30,0.31,0.30,0.36,0.32,0.30,0.51,0.30,0.31,0.37,0.30,0.47,0.29,0.28,0.40,0.30,0.31,0.39,0.49,0.30,0.00,0.00,0.00,0.00,0.01,0.00,0.09,0.01,0.19,0.00),ncol=3)
## End(Don't show)
# compare the estimated admixture proportions from 
# two different conStruct runs to determine which 
# layers in one run correspond to those in the other
match.layers.x.runs(admix.props1,admix.props2)




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
