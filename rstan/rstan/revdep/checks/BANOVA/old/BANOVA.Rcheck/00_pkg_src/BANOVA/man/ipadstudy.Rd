\name{ipadstudy}
\alias{ipadstudy}
\docType{data}
\title{Data for the study of relation between Conspicuous, Brand Usage, Self-Brand Connection and attitudes toward the brand}
\description{
  The study is a between-subjects experiment which has factor (conspicuousness: low vs. high) and one measured variable (self-brand connection). The goal is to show that conspicuous brand use negatively affects attitudes toward the user and the brand only for observers low in self-brand connection. One hundred fifty-four participants were exposed to a video manipulating conspicuous brand usage. Participants completed the study by answering several questions which are used to measure the dependent (attitude) and independent (self-brand connection) variables in the model.
}
\usage{data(ipadstudy)}
\format{
  This R object contains between-subjects variables: \$owner is an indicator variable. If the subject owns iPad or iPhone, then owner = 1. It is equal to 0 otherwise. \$age denotes the age of subjects. \$gender denotes the gender of subjects. gender = 1 if the subject is a female, 0 otherwise. \$conspic is an indicator variable related to conspicuousness. conspic = 1 if conspicuousness is high. \$selfbrand denotes the self-brand connection for Apple. \$id is the id of subjects. \$attitude denotes the attitudes towards the brand which is the continuous dependent variable. \$apple_dl is a seven-point scale variable which denotes the attitudes (dislike = 1,..., like = 7)

 
 \$ ipadstudy: 'data.frame':  154 obs. of  9 variables:\cr
 \ldots\$ id            : int  1 2 3 4 5 6 7 8 9 10 ...\cr
 \ldots\$ attitude      : num  3 5.33 5.67 5.33 6 ...\cr
 \ldots\$ owner         : num  0 0 0 1 1 0 1 0 1 0 ...\cr
 \ldots\$ age           : int  19 33 25 41 38 33 37 46 41 55 ...\cr
 \ldots\$ gender        : num  0 0 1 0 1 1 1 0 1 1 ...\cr
 \ldots\$ conspic       : num  0 1 0 1 1 0 0 1 0 1 ...\cr
 \ldots\$ selfbrand     : num  -2.304 1.696 -0.161 -0.447 0.267 ...\cr
 \ldots\$ apple_dl      : int  3 6 6 5 6 4 7 7 5 5 ...\cr
 
}

\references{ 

Ferraro,R., Kirmani, A. and Matherly, T., (2013) \emph{Look at Me! Look at Me! Conspicuous Brand Usage, Self-Brand Connection, and Dilution}. Journal of Marketing Research, Vol. 50, No. 4, pp. 477-488.

}

\examples{
data(ipadstudy)

# mean center covariates
ipadstudy$age <- ipadstudy$age - mean(ipadstudy$age) 
ipadstudy$owner <- ipadstudy$owner - mean(ipadstudy$owner )
ipadstudy$gender <- ipadstudy$gender - mean(ipadstudy$gender)

res <- BANOVA.Normal(attitude~1, ~owner + age + gender + selfbrand*conspic, 
ipadstudy, ipadstudy$id, burnin = 100, sample = 100, thin = 1 )
summary(res)
\donttest{
# use apple_dl as the dependent variable
res <- BANOVA.ordMultinomial(apple_dl~1, ~owner + age + gender + selfbrand*conspic, 
ipadstudy, ipadstudy$id, burnin = 100, sample = 100, thin = 2 )
summary(res)
table.predictions(res)
}
}