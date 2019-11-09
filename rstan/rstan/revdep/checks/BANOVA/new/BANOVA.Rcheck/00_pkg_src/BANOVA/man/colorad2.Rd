\name{colorad2}
\alias{colorad2}
\docType{data}
\title{Data for gist perception of advertising, study 2}
\description{
  Data from an experiment in which One hundred and forty eight subjects (71 men; age ranging from 21 to 28) were randomly assigned to one condition of a 2 (blur: normal, very high) x 2 (color: full color, grayscale, inverted) between-participants design. Participants were exposed to 25 ads for five brands in each of five categories. Ads were selected to be typical for the category, using the same procedure as in \code{\link{colorad}}. The product categories used were cars, financial services, food, skincare and fragrance. Images were flashed for 100 msec. and subjects were asked to identify whether the image was an ad or not, and if they identified it correctly as an ad, they were asked to indicate which category (out of five) was advertised. The total number of correct ad identifications and category identifications are used as dependent variables.
}
\usage{data(colorad2)}
\format{
  This R object contains between-subjects variables: \$B is a factor corresponding to the levels of blur (normal = 0, very high = 1), \$C1 and \$C2 are dummy variables denote 'grayscale' and 'inverted' levels of color. \$C is the original factor denote the color with 3 levels. \$ID is the ID of subjects. \$Y1 is the number of correct identifications of the 25 ads of each subject. \$Y2 is the number of correct identifications of the category, given the number of correct ad identifications.

 \$ colorad2:`data.frame':   148 obs. of  7 variables:\cr
  \ldots\$ ID  : int  1 2 3 4 5 6 7 8 9 10 ... \cr
  \ldots\$ C1  : int  0 1 1 0 0 0 0 0 1 1 ...\cr
  \ldots\$ C2  : int  0 0 0 1 1 0 0 0 0 0 ...\cr
  \ldots\$ B   : Factor w/ 2 levels "0","1": 1 1 0 0 1 0 0 1 0 1 ... \cr
  \ldots\$ Y1  : int  14 6 23 21 8 23 24 5 23 6 ... \cr
  \ldots\$ Y2  : int  2 3 8 8 2 15 10 1 13 0 ... \cr
  \ldots\$ C   : Factor w/ 3 levels "1","2","3": 1 2 2 3 3 1 1 1 2 2 ... \cr
 
}

\references{ 

Wedel, M and R. Pieters (2015). 
\emph{The Buffer Effect: The Role of Color when Advertising Exposures are Brief and Blurred}, Marketing Science, Vol. 34, No. 1, pp. 134-143.

}

\examples{ 

data(colorad2)
# factor C is effect coded
\donttest{
model_bin <- BANOVA.model('Binomial')
stanmodel_bin <- BANOVA.build(model_bin)
res0 <- BANOVA.run(Y1 ~ 1, ~ C + B + C*B, fit = stanmodel_bin, 
data = colorad2, id = colorad2$id, num_trials = as.integer(25), 
iter = 1000, thin = 1, chains = 4)
res0
# or use BANOVA.Binomial
res1 <- BANOVA.Binomial(Y1 ~ 1, ~ C + B + C*B, colorad2, colorad2$id, as.integer(25), 
burnin = 1000, sample = 1000, thin = 1) 
# or use dummy variables 
res1 <- BANOVA.Binomial(Y1 ~ 1, ~ C1 + B + C2 + C1*B + C2*B, 
colorad2, colorad2$ID, as.integer(25))
res1
summary(res1)
table.predictions(res1)

# second model
res2 <- BANOVA.Binomial(Y2 ~ 1, ~ C + B + C*B, 
colorad2, colorad2$id, colorad2$Y1)
summary(res2)
table.means(res2)
}
}