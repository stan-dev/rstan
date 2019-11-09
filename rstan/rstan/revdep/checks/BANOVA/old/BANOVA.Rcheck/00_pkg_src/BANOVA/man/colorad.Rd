\name{colorad}
\alias{colorad}
\docType{data}
\title{Data for gist perception of advertising, study 1}
\description{
  Data from an experiment in which one hundred and sixteen subjects (53 men; mean age 23, ranging from 21 to 28) were randomly assigned to one condition of a 5 (blur: normal, low, medium, high, very high) x 2 (color: full color, grayscale) between-participants, x 2 (image: typical ads, atypical ads) within-participants, mixed design. Participants were exposed to 40 images, 32 full-page ads and 8 editorial pages. There were 8 ads per product category, with 4 typical and 4 atypical ones, the categories being car, financial services, food, and skincare. Subjects were asked to identify each image being flashed for 100msec. as being an ad or not. The total number of correct ad identifications, for typical and atypical ads, are used as a dependent variable.
}
\usage{data(colorad)}
\format{
  This R object contains within-subject variable \$typic which is a factor with 2 levels "0" (typical ads) and "1"(atypical ads); between-subjects variables: \$blur which is a numerical vairable denotes 5 different levels of blur (which must be mean centered), \$blurfac is a categorical data corresponding to the levels of \$blur, \$color which is a factor with 2 levels "0"(full color) and "1"(grayscale). \$id is the ID of subjects. \$y is the number of correct identifications of the 16 ads of each subject for each level of \$typic. 

 \$ colorad:`data.frame':  474 obs. of  8 variables:\cr
  \ldots\$ id      : int  1 1 2 2 3 3 4 4 5 5 ... \cr
  \ldots\$ typic   : Factor w/ 2 levels "0","1": 0 1 0 1 0 1 0 1 0 1 ...\cr
  \ldots\$ y       : int  8 6 12 6 11 9 9 11 14 14 ... \cr
  \ldots\$ blurfac    : Factor w/ 5 levels "1","2","3","4",..: 2 2 4 4 2 2 3 3 1 1 ... \cr
  \ldots\$ color   : Factor w/ 2 levels "0","1": 1 1 0 0 0 0 0 0 1 1 ... \cr
  \ldots\$ blur: num  3.69 3.69 4.79 4.79 3.69 ... \cr
 
}

\references{ 

Wedel, M and R. Pieters (2015). 
\emph{The Buffer Effect: The Role of Color when Advertising Exposures are Brief and Blurred}, Marketing Science, Vol. 34, No. 1, pp. 134-143.

}

\examples{ 
\donttest{
data(colorad)
# Build the model
model_bin <- BANOVA.model('Binomial')
stanmodel_bin <- BANOVA.build(model_bin)
out0 = BANOVA.run(y ~ typic, ~ color*blurfac, fit = stanmodel_bin,
                   data = colorad, id = 'id', num_trials = as.integer(16), 
                   iter = 10000, thin = 10, chains = 1)
summary(out0)
# planned comparison
out0_contra = BANOVA.run(y ~ typic, ~ color*blurfac, fit = stanmodel_bin,
                          data = colorad, id = 'id', num_trials = as.integer(16), 
                          iter = 10000, thin = 10, chains = 1, 
                          contrast = list(typic = c(-1,1), 
                          blurfac = c(1,1,0,-1,-1)))
summary(out0_contra)
# floodlight analysis
set.seed(200)
out0_f = BANOVA.run(y ~ typic, ~ color*blur, fit = stanmodel_bin,
                     data = colorad, id = 'id', num_trials = as.integer(16), 
                     iter = 10000, thin = 10, chains = 1)
BANOVA.floodlight(out0_f, var_factor = 'typic', var_numeric = 'blur')
}

}