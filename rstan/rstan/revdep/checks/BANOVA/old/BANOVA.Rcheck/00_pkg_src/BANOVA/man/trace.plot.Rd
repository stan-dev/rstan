\name{trace.plot}
\alias{trace.plot}
\title{Function to plot the trace of parameters}
\description{
Function to plot the trace of all coefficients/parameters. The plots can be saved as a pdf file.
}
\usage{
trace.plot(x, save = FALSE)
}
\arguments{
  \item{x}{the object from BANOVA.*}

  \item{save}{whether to save the trace plot as a pdf file, the default is FALSE}
}

\examples{
data(goalstudy)
res1 <- BANOVA.Normal(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 1000, sample = 1000, thin = 2) 
\donttest{
# or use BANOVA.run
res1 <- BANOVA.run(bid~1, ~progress*prodvar, model_name = "Normal", 
data = goalstudy, id = goalstudy$id, iter = 1000, thin = 1, chains = 2) 
}
trace.plot(res1)

}

