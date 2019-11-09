\name{condstudy}
\alias{condstudy}
\docType{data}
\title{Data for the study of how brand attitudes were influenced by showing brands together with pleasant pictures}
\description{
  The study investigated how brand attitudes were influenced by showing brands together with pleasant pictures. Attitude change via conditioning can result from either a direct transfer of affect from the picture to the brand, or from an indirect association of the brand and the picture in memory. In Sweldens' et al. (2010) experiment 1, indirect conditioning was implemented by presenting a brand repeatedly with the same picture, direct conditioning by presenting it simultaneously with different pictures. The pictures used were either neutral or positive. This study involved a mixed design, with a within-subject factor (cond = neutral, positive), and a between-subject factor (type = indirect, direct), as well as a within-subject mediator. Although the original mediation hypotheses are more intricate, here the mediation of the conditioning effect is investigated by measurements of attitudes towards the pictures that were shown with the brands (pict).
}
\usage{data(condstudy)}
\format{
  This R object contains a between-subjects variable: type, which denotes a between-subject moderator. It has two levels, "indirect" and "direct". In the "indirect" condition the brands were shown with the same images, in the indirect condition the brands were shown with different images; Within-subject variables: cond, a within-subject factor with 2 levels: "pos", and "xneu", which indicates whether each brand was shown with a neutral (xneu) or positive (pos) emotional image. pict,  a within-subject mediator variable measuring the valence (positive/negative) of the emotional image the respondent remembers the brand to have been shown with. att, a dependent variable which denotes the ratings of attitudes toward brands. 

 
 \$ condstudy: 'data.frame':  888 obs. of  5 variables:\cr
 \ldots\$ id  : int  2 2 2 2 2 2 3 3 3 3 ... \cr
 \ldots\$ att : num  2.94 2.44 3.44 1.67 1.67 ... \cr
 \ldots\$ cond: Factor w/ 2 levels "pos","xneu": 1 1 1 2 2 2 1 1 1 2 ... \cr
 \ldots\$ type: Factor w/ 2 levels "direct","indirect": 2 2 2 2 2 2 2 2 2 2 ... \cr
 \ldots\$ pict: int  6 7 6 2 4 5 9 3 2 5 ... \cr
 
}

\references{ 

Sweldens, S., Osselaer, S. and Janiszewski, C. (2010) \emph{Evaluative Conditioning Procedures and the Resilience of Conditioned Brand Attitudes}. Journal of Consumer Research,  Vol. 37.

Wedel, M. and Dong, C. (2016) \emph{BANOVA: Bayesian Analysis of Variance for Consumer Research}. Submitted.

}

\examples{
data(condstudy)
\donttest{
model <- BANOVA.model('Normal')
stanmodel <- BANOVA.build(model)
out2 <- BANOVA.run(att~cond+pict, ~type, fit = stanmodel, data = condstudy, 
                  id = condstudy$id, iter = 500, thin = 1, chains = 2)    
conv.diag(out2)
summary(out2) 
table.predictions(out2)  
BANOVA.floodlight(out2, var_factor = 'type', var_numeric = 'pict')
cat(out2$model_code)

out3 <- BANOVA.run(pict~cond, ~type, fit = stanmodel, data = condstudy, 
                  id = condstudy$id, iter = 500, thin = 1, chains = 2)
conv.diag(out3)
summary(out3)
BANOVA.mediation(out2, out3, xvar='cond', mediator='pict')
}
}