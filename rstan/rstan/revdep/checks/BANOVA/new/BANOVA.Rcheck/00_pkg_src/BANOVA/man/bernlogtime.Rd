\name{bernlogtime}
\alias{bernlogtime}
\docType{data}
\title{Data for analysis of effects of typicality, blur and color on gist perception of ads}
\description{
  Data from a mixed design experiment, where respondents were exposed to 32 ads, for 100 millisec. The ads were either typical or atypical (typical: 1 or 2). Respondents were exposed to ads that were either in full color or black-and-white (color: 1 or 2), and at different levels of blur (1=normal,5 = very high blur). These are between-subjects factors. The dependent variables are the response 0/1, and the response time. Typicality is a within-subjects variable.
}
\usage{data(bernlogtime)}
\format{
  This R object contains within-subject variable: \$typical is a factor with 2 levels "0" (typical ads) and "1"(atypical ads); between-subjects variables: \$blur is a factor with two levels (1=normal,5 = very high blur). \$color denotes a factor with 2 levels "1"(full color) and "2"(grayscale). \$subject is the ID of subjects. \$response denotes if the ad is correctly identified. \$logtime is the response time.

 \$bernlogtime: 'data.frame':  3072 obs. of  6 variables:\cr
 \ldots\$ subject : int  5 5 5 5 5 5 5 5 5 5 ...\cr
 \ldots\$ typical : Factor w/ 2 levels "1","2": 1 2 1 1 1 2 2 2 2 1 ...\cr
 \ldots\$ blur    : Factor w/ 2 levels "1","5": 1 1 1 1 1 1 1 1 1 1 ...\cr
 \ldots\$ color   : Factor w/ 2 levels "1","2": 2 2 2 2 2 2 2 2 2 2 ...\cr
 \ldots\$ response: int  1 1 1 1 1 1 1 1 1 1 ...\cr
 \ldots\$ logtime : num  0.977 1.73 1.784 1 1.149 ...\cr
 
}

\references{ 

Wedel, M and R. Pieters (2015). 
\emph{The Buffer Effect: The Role of Color when Advertising Exposures are Brief and Blurred}, Marketing Science, Vol. 34, No. 1, pp. 134-143.

}

\examples{
data(bernlogtime)
\donttest{
# model using the dependent variable : log of the response time(logtime) 
res1 <- BANOVA.Normal(logtime~typical, ~blur + color, bernlogtime, 
bernlogtime$subject, burnin = 1000, sample = 1000, thin = 1)
summary(res1)
table.predictions(res1)

# model using the dependent variable : response
res2 <- BANOVA.Bernoulli(response~typical, ~blur + color, bernlogtime, 
bernlogtime$subject, burnin = 1000, sample = 1000, thin = 1)
summary(res2)
table.predictions(res2)
}
}