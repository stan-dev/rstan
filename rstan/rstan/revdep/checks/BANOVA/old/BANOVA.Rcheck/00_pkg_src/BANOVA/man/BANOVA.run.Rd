\name{BANOVA.run}
\alias{BANOVA.run}
\alias{predict.BANOVA}
\alias{print.BANOVA}
\alias{summary.BANOVA}
\title{Estimation of BANOVA models}
\description{
\code{BANOVA.run} implements Hierarchical Bayesian ANOVA models using 'Stan'
}
\usage{
BANOVA.run(l1_formula = "NA", l2_formula = "NA", fit = NULL, model_name = 'NA', 
dataX = NULL, dataZ = NULL, data = NULL, y_value = NULL, id, iter = 2000, 
num_trials = 1, contrast = NULL, ...)
\method{summary}{BANOVA}(object, ...)
\method{predict}{BANOVA}(object, newdata = NULL, Xsamples = NULL, Zsamples = 
                        NULL, ...)
\method{print}{BANOVA}(x, ...)
}
\arguments{
  \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included. If NULL, the single level model is used}

  \item{fit}{a fitted BANOVA models, an object of class \code{"BANOVA.build"}, default NULL which needs compilation}
  
  \item{model_name}{a character string in c('Normal', 'T', 'Bernoulli', 'Binomial', 'Poisson', 'ordMultinomial', 'Multinomial')}

  \item{dataX}{a list of data frames(each corresponds to the choice set of each observation) that includes all covariates and factors, for the Multinomial model only, default NULL}
  
  \item{dataZ}{a data frame(long format) that includes all level 2 covariates and factors, for the Multinomial model only, default NULL}
  
  \item{data}{a data.frame in long format including all features in level 1 and level 2(covariates and categorical factors) and responses, default NULL}

  \item{id}{subject ID of each response unit}
  
  \item{y_value}{choice responses, 1,2,3..., for the Multinomial model only, default NULL}

  \item{iter}{target samples in the 'Stan' algorithm after thinning, default 2000}
  
  \item{num_trials}{the number of trials of each
  observation(=1, if it is bernoulli), the type is forced to be 'integer', for the Binomial model only, default 0}
  
  \item{contrast}{a list of contrasts for planned comparisons, default: effect coding (NULL value)}

  \item{object}{an object of class \code{BANOVA} (returned by \code{BANOVA.run})}
  
  \item{x}{an object of class \code{BANOVA} (returned by \code{BANOVA.run})}

  \item{newdata}{test data, either a matrix, vector or a
  data frame. It must have the same format with the original data (the same column number)}
  \item{Xsamples}{a list of sample data frames(each corresponds to the choice set of each observation) that includes all covariates and factors, for the Multinomial model only, default NULL}
  \item{Zsamples}{a data frame(long format) that includes all level 2 covariates and factors, for the Multinomial model only, default NULL}
  \item{\dots}{additional arguments, for \code{BANOVA.run}, it can include standard 'Stan' arguments, e.g. warmup, thin, chains, etc., see \link[rstan]{sampling} for more details, for other functions, ignored currently}
}

\value{
  \code{BANOVA.run} returns an object of class \code{"BANOVA"}. The returned object is a list containing:
  
  \item{anova.table}{table of effect sizes \code{\link{BAnova}}}
  \item{coef.tables}{table of estimated coefficients}
  \item{pvalue.table}{table of p-values \code{\link{table.pvalues}}}
  \item{dMatrice}{design matrices at level 1 and level 2}
  \item{samples_l1_param}{posterior samples of level 1 parameters}
  \item{samples_l2_param}{posterior samples of level 2 parameters}
  \item{samples_l2_sigma_param}{posterior samples of level 2 standard deviations}
  \item{samples_cutp_param}{posterior samples of cutpoints}
  \item{data}{original data.frame}
  \item{mf1}{model.frame of level 1}
  \item{mf2}{model.frame of level 2}
  \item{model_code}{'Stan' code}
  \item{single_level}{if this is a single level model}
  \item{stan_fit}{fitted samples}
  \item{model_name}{the name of the model}
  \item{contrast}{contrasts for planned comparisons}
  \item{new_id}{id values coded in 1,2,3,...}
  \item{old_id}{original id values}
}

\examples{
# Use the ipadstudy data set
data(ipadstudy)
\donttest{
# build the BANOVA model first so that it can be reused
model <- BANOVA.model('Normal', single_level = FALSE)
banova_model <- BANOVA.build(model)
res_1 <- BANOVA.run(attitude~1, ~owner + age + gender + selfbrand*conspic, 
fit = banova_model, data = ipadstudy, id = ipadstudy$id, iter = 2000, 
thin = 5, chains = 2)
summary(res_1)
res_2 <- BANOVA.run(attitude~1, ~selfbrand*conspic, 
fit = banova_model, data = ipadstudy, id = ipadstudy$id, iter = 2000, 
thin = 5, chains = 2)

# or call the function directly without specifying the fit argument
# but it needs compilation
res_1 <- BANOVA.run(attitude~1, ~owner + age + gender + selfbrand*conspic,
model_name = 'Normal', data = ipadstudy, id = ipadstudy$id, iter = 2000, 
thin = 5, chains = 2)
}

}


