#!/usr/bin/r -t
#
# Copyright (C) 2011 - 2015  Douglas Bates, Dirk Eddelbuettel and Romain Francois
#
# This file is part of RcppEigen
#
# RcppEigen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppEigen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

.setUp <- function(){
    suppressMessages(require(datasets))
    suppressMessages(require(RcppEigen))
}

test.fastLm <- function() {
    data(trees, package="datasets")
    flm0 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 0L,
                  PACKAGE="RcppEigen")
    flm1 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 1L,
                  PACKAGE="RcppEigen")
    flm2 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 2L,
                  PACKAGE="RcppEigen")
    flm3 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 3L,
                  PACKAGE="RcppEigen")
    flm4 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 4L,
                  PACKAGE="RcppEigen")
    flm5 <- .Call("RcppEigen_fastLm_Impl",
                  cbind(1, log(trees$Girth)),
                  log(trees$Volume), 5L,
                  PACKAGE="RcppEigen")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)
    fitCoef <- unname(coef(fit))
    fitStdErr <- unname(coef(summary(fit))[, "Std. Error", drop = TRUE])
    checkEquals(flm0$coefficients, fitCoef, msg="fastLm0.coef")
    checkEquals(flm0$se, fitStdErr, msg="fastLm0.stderr")
    checkEquals(flm1$coefficients, fitCoef, msg="fastLm1.coef")
    checkEquals(flm1$se, fitStdErr, msg="fastLm1.stderr")
    checkEquals(flm2$coefficients, fitCoef, msg="fastLm2.coef")
    checkEquals(flm2$se, fitStdErr, msg="fastLm2.stderr")
    checkEquals(flm3$coefficients, fitCoef, msg="fastLm3.coef")
    checkEquals(flm3$se, fitStdErr, msg="fastLm3.stderr")
    checkEquals(flm4$coefficients, fitCoef, msg="fastLm0.coef")
    checkEquals(flm4$se, fitStdErr, msg="fastLm0.stderr")
    checkEquals(flm5$coefficients, fitCoef, msg="fastLm0.coef")
    checkEquals(flm5$se, fitStdErr, msg="fastLm0.stderr")
}


test.fastLm.formula <- function() {
    data(trees, package="datasets")
    flm <- fastLm(log(Volume) ~ log(Girth), data=trees)
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(flm$coefficients, coef(fit), msg="fastLm.formula.coef")
    checkEquals(as.numeric(flm$se), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.formula.stderr")
}

