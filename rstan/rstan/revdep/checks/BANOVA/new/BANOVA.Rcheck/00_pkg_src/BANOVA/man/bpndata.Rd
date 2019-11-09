\name{bpndata}
\alias{bpndata}
\docType{data}
\title{Eye-movement data for analysis of print ad designs}
\description{
Data were collected in an experimental study in which 88 participants freely paged through a magazine at home or in a waiting room. While flipping through pages at their own pace, participants' eye-movements were recorded with infra-red corneal reflection eye-tracking methodology. In a subsequent memory task, participants were asked to identify the target brand in the ad as soon as possible by touching the correct brand name on the screen. \bold{Accuracy} (accurate=1, inaccurate =0) of brand memory and \bold{response time} were recorded for each ad and participant. 

}
\usage{data(bpndata)}
\format{
  This R object contains 3080 observations in the data (35 ads x 88 participants). The goal is to examine the effects of several ad design variables on both eye movements and memory. The variables include:

1.  RESPONDENT_ID: ID number of a respondent;\cr
2.  AD_ID: ID number of an ad;\cr
3.	PAGE_NUM: page number in the magazine where an ad appears (1,2,3,...);\cr
4.	PAGE_POS: the right-side vs. left-side position on a page, 1 = right, 0 = left;\cr
5.	PIC_FIX: fixation count of the pictorial element (0, 1, 2, 3, ...);\cr
6.	PIC_SIZE: surface size of the pictorial element, in inches2;\cr
7.	RECALL_ACCU: whether a respondent accurately recalls the brand name, 1= yes, 0 = no;\cr
8.	RECALL_TIME: the time it takes a respondent to answer the brand recall question, in seconds.\cr

 \$ bpndata: 'data.frame':  3080 obs. of  8 variables:\cr
 \ldots\$ RESPONDENT_ID: int  1 1 1 1 1 1 1 1 1 1 ...\cr
 \ldots\$ AD_ID        : int  1 2 3 4 5 6 7 8 9 10 ...\cr
 \ldots\$ PAGE_NUM     : int  2 5 6 11 13 14 17 18 21 22 ...\cr
 \ldots\$ PAGE_POS     : int  0 1 0 1 1 0 1 0 1 0 ...\cr
 \ldots\$ PIC_FIX      : int  0 2 1 1 1 2 0 3 3 8 ...\cr
 \ldots\$ PIC_SIZE     : num  74.2 52.6 77.6 71.4 52.4 ...\cr
 \ldots\$ RECALL_ACCU  : int  0 0 0 0 0 0 1 1 0 0 ...\cr
 \ldots\$ RECALL_TIME  : num  2.56 1.04 2.76 2.8 2.28 2.32 2.04 2.04 2.48 0.6 ...\cr
}

\references{ 

Wedel, M. and Pieters, R. (Autumn, 2000). \emph{Eye Fixations on Advertisements and Memory for Brands: A Model and Findings}, Marketing Science, Vol. 19, No. 4, pp. 297-312

}

\examples{ 
data(bpndata)
# within-subjects model using the dependent variable : PIC_FIX
\donttest{
model <- BANOVA.model('Poisson')
stanmodel <- BANOVA.build(model)
res0 <- BANOVA.run(PIC_FIX ~ PIC_SIZE + PAGE_NUM + PAGE_POS, ~1, 
fit = stanmodel, data = bpndata, id = bpndata$RESPONDENT_ID, 
iter = 200, thin = 1, chains = 2)
res0
# or
res1 <- BANOVA.Poisson(PIC_FIX ~ PIC_SIZE + PAGE_NUM 
+ PAGE_POS, ~1, bpndata, bpndata$RESPONDENT_ID, burnin = 1000, sample = 1000, thin = 1)
res1

# within-subjects model using the dependent variable : RECALL_ACCU
model_bern <- BANOVA.model('Bernoulli')
stanmodel_bern <- BANOVA.build(model_bern)
res2 <- BANOVA.run(RECALL_ACCU ~ RECALL_TIME + PAGE_NUM + PAGE_POS, ~1, 
fit = stanmodel_bern, data = bpndata, id = bpndata$RESPONDENT_ID, 
iter = 200, thin = 1, chains = 2)
res2
# or
res3 <- BANOVA.Bernoulli(RECALL_ACCU ~ RECALL_TIME + PAGE_NUM 
+ PAGE_POS, ~1, bpndata, bpndata$RESPONDENT_ID, burnin = 1000, sample = 1000, thin = 1)
res3
}
}