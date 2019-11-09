## ---- cache=FALSE, include=FALSE-----------------------------------------
library(knitr)
library(rpf)
library(ggplot2)
library(reshape2)
library(gridExtra)
opts_chunk$set(echo=TRUE)

## ------------------------------------------------------------------------
spec <- list()
spec[1:6] <- rpf.grm(factors=2)
gen.param <- sapply(spec, rpf.rparam)
colnames(gen.param) <- paste("i", 1:ncol(gen.param), sep="")
gen.param[2,] <- c(0,0,.5,.5,1,1)

resp <- rpf.sample(1000, spec, gen.param)

# hide latent factor that we don't know about
tspec <- list()
tspec[1:length(spec)] <- rpf.grm(factors=1)

grp <- list(spec=tspec, param=gen.param[-2,], mean=c(0), cov=diag(1), data=resp)

ChenThissen1997(grp)

## ------------------------------------------------------------------------
(got <- SitemFit(grp))

## ---- echo=FALSE,fig.height=2.5------------------------------------------
SSplot <- function(sout, itemName, showSampleSize=TRUE) {
    s1 <- sout[[itemName]]
    obs <- s1$orig.observed
    ex <- s1$orig.expected
    rowTotal <- apply(obs, 1, sum)
    mask <- rowTotal > 0
    obs <- (obs / rowTotal)[mask,]
    ex <- (ex / rowTotal)[mask,]
    ss <- data.frame(sscore=as.numeric(names(rowTotal)), n=rowTotal)
    both <- rbind(cbind(type="expected", melt(ex)),
                  cbind(type="observed", melt(obs)))
    both$outcome <- factor(both$outcome, colnames(obs))
    plot <- ggplot(both, aes(x=sumScore, y=value)) + facet_wrap(~type) + ylim(0,1) +
        labs(y="probability", title=itemName)
    guide.style <- guide_legend(keywidth=.1, keyheight=.5, direction = "horizontal", title.position = "top",
                                label.position="bottom", label.hjust = 0.5, label.vjust = .5,
                                label.theme = element_text(angle = 90, size=8))
    plot <- plot + geom_line(aes(color=outcome)) + guides(color = guide.style)
    if (showSampleSize) {
        plot <- plot + geom_text(data=ss, aes(label=n, x=sscore), y = 1, size=2, angle=90)
    }
    plot
}
SSplot(got, "i1")
SSplot(got, "i2")
SSplot(got, "i3")
SSplot(got, "i4")
SSplot(got, "i5")
SSplot(got, "i6")

## ------------------------------------------------------------------------
(got <- sumScoreEAP(grp))

## ---- echo=FALSE,fig.height=2.5------------------------------------------
got <- sumScoreEAPTest(grp)
df <- data.frame(score=as.numeric(rownames(got$tbl)),
            expected=got$tbl[,'p'] * got$n, observed=got$observed)
df <- melt(df, id="score", variable.name="source", value.name="n")
ggplot(df, aes(x=score, y=n, color=source)) + geom_line()

## ------------------------------------------------------------------------
  data(science)
  spec <- list()
  spec[1:25] <- rpf.nrm(outcomes=3, T.c = lower.tri(diag(2),TRUE) * -1)
  param <- rbind(a=1, alf1=1, alf2=0,
        gam1=sfif$MEASURE + sfsf[sfsf$CATEGORY==1,"Rasch.Andrich.threshold.MEASURE"],
        gam2=sfif$MEASURE + sfsf[sfsf$CATEGORY==2,"Rasch.Andrich.threshold.MEASURE"])
  colnames(param) <- sfif$NAME
  iorder <- match(sfif$NAME, colnames(sfpf))
  responses <- sfpf[,iorder]
  rownames(responses) <- sfpf$NAME
  
  rpf.1dim.fit(spec, param, responses, sfpf$MEASURE, 2, wh.exact=TRUE)
  head(rpf.1dim.fit(spec, param, responses, sfpf$MEASURE, 1, wh.exact=TRUE))

