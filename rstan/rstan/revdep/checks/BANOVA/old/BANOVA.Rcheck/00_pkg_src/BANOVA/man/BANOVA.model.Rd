\name{BANOVA.model}
\alias{BANOVA.model}
\title{Extract BANOVA models}
\description{
\code{BANOVA.model} extracts BANOVA models from the package.
}
\usage{
BANOVA.model(model_name, single_level = F)
}
\arguments{
  \item{model_name}{a character string in c('Normal', 'T', 'Bernoulli', 'Binomial', 'Poisson', 'ordMultinomial', 'Multinomial')}

  \item{single_level}{if the model is a single level model, default False}
}

\value{
  \code{BANOVA.model} returns an object of class \code{"BANOVA.model"}. The returned object is a list containing:
  
  \item{model_code}{the model code of the extracted model}
  \item{model_name}{the model name}
  \item{single_level}{if the model is a single level model}

}

\examples{
\donttest{
model <- BANOVA.model('Poisson', single_level = FALSE) 
cat(model$model_code)
}
}



