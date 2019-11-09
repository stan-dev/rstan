\name{choicedata}
\alias{choicedata}
\docType{data}
\title{Household Panel Data on Margarine Purchases}
\description{
  Panel data on purchases of margarine by 204 households. 
  Demographic variables are included.
}
\usage{data(choicedata)}
\format{
  This is an R object that contains within-subjects variables and between-subjects variables:

 \$ choicePrice:`data.frame':  1500 obs. of  13 variables:\cr
  \ldots \$ hhid    : int  2100016 2100016 2100016 2100016 \cr
  \ldots \$ choice  : int  1 1 1 1 1 4 1 1 4 1 \cr
  
  Within-subject variables: \cr
  
  \ldots \$ PPk\_Stk : num 0.66 0.63 0.29 0.62 0.5 0.58 0.29 ...\cr
  \ldots \$ PBB\_Stk : num 0.67 0.67 0.5 0.61 0.58 0.45 0.51 ...\cr
  \ldots \$ PFl\_Stk : num 1.09 0.99 0.99 0.99 0.99 0.99 0.99 ...\cr
  \ldots \$ PHse\_Stk: num 0.57 0.57 0.57 0.57 0.45 0.45 0.29 ...\cr
  \ldots \$ PGen\_Stk: num 0.36 0.36 0.36 0.36 0.33 0.33 0.33 ...\cr
  \ldots \$ PSS\_Tub : num 0.85 0.85 0.79 0.85 0.85 0.85 0.85 ...\cr

  Pk is Parkay; BB is BlueBonnett, Fl is Fleischmanns, Hse is house,
  Gen is generic, SS is Shed Spread.  \_Stk indicates 
  stick, \_Tub indicates Tub form.
  
  Between-subject variables: \cr
  
  \ldots \$ Income   : num 32.5 17.5 37.5 17.5 87.5 12.5 ... \cr
  \ldots \$ Fam\_Size : int 2 3 2 1 1 2 2 2 5 2 ...\cr
  \ldots \$ college  : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...\cr
  \ldots \$ whtcollar: Factor w/ 2 levels "0","1": 0 0 0 0 0 0 0 1 1 1 ...\cr
  \ldots \$ retired  : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...\cr

}
\details{
  choice is a multinomial indicator of one of the 6 brands (in order listed under format). 
  All prices are in $.
}
\source{
  Allenby, G. and Rossi, P. (1991), \emph{Quality Perceptions and Asymmetric Switching Between Brands}, \emph{Marketing Science}, Vol. 10, No.3, pp. 185-205.
}
\references{
  Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi et al.\cr
  \url{http://www.perossi.org/home/bsm-1}
}
\examples{ 
data(choicedata)
# generate dataX(convert the within-subjects variables to a list)
dataX <- list()
for (i in 1:nrow(choicedata)){
  logP <- as.numeric(log(choicedata[i,3:8]))
  # Note: Before the model initialization, all numeric variables(covariates) 
  # must be mean centered
  dataX[[i]] <- as.data.frame(logP) - mean(logP)
}
dataZ <- choicedata[,9:13]
\donttest{
res <- BANOVA.Multinomial(~ logP, ~ college, dataX, dataZ, choicedata$choice, 
choicedata$hhid, burnin = 1000, sample = 1000, thin = 1)
summary(res)
predict(res,dataX[1:4], dataZ[1:4,])
}
}