\name{BAnova}
\alias{BAnova}
\title{Function to print the table of effect sizes}
\description{
The analysis of variance is performed at level 1 (for the single level model) and level 2 equation of the Bayesian ANOVA see \code{\link{BANOVA-package}}. This makes it possible to capture the effects of level-1 or level-2 variables on the heterogeneity distribution of subjects, and compute sums of squares and effect sizes.

}
\usage{
BAnova(x)
}
\arguments{
  \item{x}{the object from BANOVA.*}
}

\details{
Measures of effect size in regression are measures of the degree of association between an effect (e.g., a main effect, an interaction, a linear contrast) and the dependent variable. They can be considered as the correlation between a categorical factor(effect) and the dependent variable. They are usually interpreted as the proportion of variance in the dependent variable that is attributable to each effect. In the package, partial Eta squared is calculated and displayed. It is defined as follows, \cr

\eqn{\eta^2= \frac{(SS_{effect})}{(SS_{effect} + SS_{error})}} \cr

Where:
SS_effect= the sums of squares for the effect of interest \cr
SS_error= the sums of squares for the error in the regression.

This equation is evaluated at each draw of the parameters, which allows for the calculation of not only the posterior mean, but also the credible interval of the effect size. 

}

\references{ 
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Lakens, D. (2013) 
  \emph{Calculating and Reporting Effect Sizes to Facilitate Cumulative Science: A Practical Primer for T-tests and ANOVAs}, Frontiers in Psychology, Vol. 4, pp.863.
  
  Gelman, A. and  Pardoe, I. (2006)
  \emph{Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models}, TECHNOMETRICS, Vol. 48, NO. 2.
  
}

\examples{

data(goalstudy)
res1 <- BANOVA.Normal(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 1000, sample = 1000, thin = 2) 
BAnova(res1)
\donttest{
# or use BANOVA.run based on 'Stan'
res2 <- BANOVA.run(bid~1, ~progress*prodvar,  model_name = 'Normal', 
data = goalstudy, id = goalstudy$id, iter = 1000, chains = 2)
BAnova(res2)
}

}

