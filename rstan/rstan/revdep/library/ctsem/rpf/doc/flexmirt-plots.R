## ---- cache=FALSE, include=FALSE-----------------------------------------
library(knitr)
library(rpf)
library(ggplot2)
library(reshape2)
library(gridExtra)
opts_chunk$set(echo=FALSE)

## ----echo=TRUE-----------------------------------------------------------
small <- structure(list(param = structure(c(1, 1, 0, 0, -0.5789195, -2.412259,  -1.3471789, 1, 1, 0, 0, 1.0983234, -2.0991327, -2.9482965, 1,  1, 0, 0, 0.4078264, -0.9824549, -1.5905594, 1, 1, 0, 0, -1.0650001,  -0.2100243, -3.2034577), .Dim = c(7L, 4L), .Dimnames = list(NULL,      c("X2", "X6", "X7", "X10"))), mean = 0, cov = structure(9.8238066, .Dim = c(1L,  1L))), .Names = c("param", "mean", "cov"))
small$spec <- list()
small$spec[1:4] <- rpf.nrm(outcomes=4, T.c= lower.tri(diag(3),TRUE) * -1)

## ----fig.height=5--------------------------------------------------------
width <- 5
small$icc <- list()
small$iif <- list()

tcc <- expand.grid(theta=seq(-width,width,.1), score=0)
tic <- expand.grid(theta=seq(-width,width,.1), info=0)

for (ix in 1:length(small$spec)) {
    name <- colnames(small$param)[ix]
    ii <- small$spec[[ix]]
    ii.p <- small$param[,ix]
    grid <- expand.grid(theta=seq(-width,width,.1))
    grid <- cbind(grid, t(rpf.prob(ii, ii.p, grid$theta)))
    tcc$score <- tcc$score + c(0:(ii@outcomes-1) %*% rpf.prob(ii, ii.p, grid$theta))
    colnames(grid) <- c("theta", paste0("k", 0:3))
    grid2 <- melt(grid, id.vars=c("theta"), variable.name="category", value.name="p")

    small$icc[[ix]] <- ggplot(grid2, aes(theta, p, color=category)) + geom_line() +
            ggtitle(paste("Item", name)) + ylim(0,1) + xlim(-width, width)

    grid <- expand.grid(theta=seq(-width,width,.1))
    grid$info <- rpf.info(ii, ii.p, t(grid$theta))
    tic$info <- tic$info + grid$info
    small$iif[[ix]] <- ggplot(grid, aes(theta, info)) + geom_line() +
            ggtitle(paste("Item", name)) + xlim(-width, width)
  }
do.call(grid.arrange, c(small$icc, ncol=2))

## ----fig.height=5--------------------------------------------------------
do.call(grid.arrange, c(small$iif, ncol=2))

## ----fig.height=2.5------------------------------------------------------
tcc.plot <- ggplot(tcc, aes(theta, score)) + geom_line() +
  ggtitle("Test Characteristic Curve") + xlim(-width, width)
tic.plot <- ggplot(tic, aes(theta, info)) + geom_line() +
  ggtitle("Test Information Curve") + xlim(-width, width)
grid.arrange(tcc.plot, tic.plot, ncol=2)

