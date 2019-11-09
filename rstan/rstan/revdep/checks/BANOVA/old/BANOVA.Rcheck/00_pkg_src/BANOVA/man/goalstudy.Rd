\name{goalstudy}
\alias{goalstudy}
\docType{data}
\title{Data for the study of the impact of the variety among means on motivation to pursue a goal}
\description{
  The study investigated how the perceived variety (high vs. low) among products, as means to a subjects' goal, affects their motivation to pursue that goal. The hypothesis was that only when progress toward a goal is low, product variety increases motivation to pursue the goal. In the study,  one hundred and five subjects were randomly assigned to conditions in a 2 (goal progress: low vs. high) by 2 (variety among means: low vs. high) between-subjects design. The final goal was a "fitness goal", and the products used were protein bars; variety was manipulated by asking subjects to think about how the products were similar (low) or different (high); goal progress was primed by asking subjects questions regarding the frequency of their recent workouts on low (0,1,...,5 or more) versus high (5 or less, 6,7,..., 10) frequency scales. Subjects were asked questions regarding the similarity of protein bars, and the bid they were willing to make for the bars, used as dependent variables in the study. 
}
\usage{data(goalstudy)}
\format{
  This R object contains between-subjects variables: progress, which denotes the progress toward a goal (1:low , 2: high ); prodvar, which denotes the amount of variety within the means to goal attainment (1:low , 2:high); perceivedsim, which is a seven-point scale dependent variable measuring the perceived similarity of the set of products (1 = not at all similar, 7 = very similar); and bid which denotes the amount that subjects would be willing to pay for the products . 

 
 \$ goalstudy: 'data.frame':  105 obs. of  5 variables:\cr
 \ldots\$ id            : int  1 2 3 4 5 6 7 8 9 10 ...\cr
 \ldots\$ perceivedsim  : int  5 7 2 2 5 5 5 4 5 7 ... \cr
 \ldots\$ progress      : Factor w/ 2 levels "1","2": 1 1 2 2 2 1 2 1 2 1 ... \cr
 \ldots\$ prodvar      : Factor w/ 2 levels "1","2": 2 1 2 1 1 1 1 2 1 1 ...\cr
 \ldots\$ bid         : num  5 0 1 15 3 10 5 4.5 3 0.75 ...\cr
 
}

\references{ 

Etkin, J. and Ratner, R. (2012) \emph{The Dynamic Impact of Variety among Means on Motivation}. Journal of Consumer Research,  Vol. 38, No. 6, pp. 1076 - 1092.

}

\examples{
data(goalstudy)

res1 <- BANOVA.Normal(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 1000, sample = 1000, thin = 1) 
\donttest{
# use BANOVA.run based on 'Stan'
res1 <- BANOVA.run(bid~1, ~progress*prodvar, model_name = "Normal", 
data = goalstudy, id = goalstudy$id, iter = 1000, thin = 1, chains = 2) 
conv.diag(res1)
# using pairs, it shows that the chains are not converged
# the level 1 and level 2 variances are not identified since the data 
# has only one obervation per subject, it is better to use single level models
pairs(res1, pars = c("beta1[1,1]", "beta2[1,1]", "tau_ySq", "tau_beta1Sq"))
# single level model
res1 <- BANOVA.run(bid~progress*prodvar, model_name = "Normal", 
data = goalstudy, id = goalstudy$id, iter = 1000, thin = 1, chains = 2) 
BAnova(res1)
table.pvalues(res1)
trace.plot(res1)
table.predictions(res1)

res1_t <- BANOVA.T(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 1000, sample = 1000, thin = 1) 

res2 <- BANOVA.ordMultinomial (perceivedsim~1, ~progress*prodvar,  goalstudy,
goalstudy$id, burnin = 1000, sample = 1000, thin = 1)

summary(res2)
table.predictions(res2)
# print 'JAGS' code
cat(res2$JAGSmodel)
}
}