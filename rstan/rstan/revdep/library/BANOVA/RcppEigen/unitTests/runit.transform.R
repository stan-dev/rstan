#
# Copyright (C) 2012 - 2013  Douglas Bates, Dirk Eddelbuettel and Romain Francois
#
# This file is part of RcppEigen.
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
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

incl <- '
typedef Eigen::ArrayXd                   Ar1;
typedef Eigen::Map<Ar1>                 MAr1;
typedef Eigen::ArrayXXd                  Ar2;
typedef Eigen::Map<Ar2>                 MAr2;
typedef Eigen::MatrixXd                  Mat;
typedef Eigen::Map<Mat>                 MMat;
typedef Eigen::VectorXd                  Vec;
typedef Eigen::Map<Vec>                 MVec;
'

definitions <- list(
    "ar1_unbounded" = list(signature(x_="numeric"),
    '
    MAr1           x(as<MAr1>(x_));

    return List::create(Named("abs",    x.abs()),
                        Named("abs2",   x.abs2()),
                        Named("exp",    x.exp()),
                        Named("cos",    x.cos()));
    '),
    "ar2_unbounded" = list(signature(X_="matrix"),
    '
    MAr2           X(as<MAr2>(X_));

    return List::create(Named("abs",    X.abs()),
                        Named("abs2",   X.abs2()),
                        Named("exp",    X.exp()),
                        Named("cos",    X.cos()));
    ')
    )

.setUp <- function() {
    suppressMessages(require(inline))
    suppressMessages(require(RcppEigen))
    cxxargs <- ifelse(Rcpp:::capabilities()[["initializer lists"]],
                      "-std=c++0x","")
    tests <- ".rcppeigen.trans"
    if( ! exists( tests, globalenv() )) {
        fun <- RcppEigen:::compile_unit_tests(definitions,
                                              includes=incl,
                                              cxxargs = cxxargs)
        names(fun) <- names(definitions)
        assign(tests, fun, globalenv())
    }
}

test.transformationAr1 <- function() {
    set.seed(1234321)
    x <- rnorm(10L)

    res <- .rcppeigen.trans$ar1_unbounded(x)
    checkEquals(res$abs,  abs(x))
    checkEquals(res$abs2, x * x)
    checkEquals(res$exp,  exp(x))
    checkEquals(res$cos,  cos(x))
}

test.transformationAr2 <- function() {
    set.seed(1234321)
    X <- matrix(rnorm(100L), nrow = 10, ncol = 10)

    res <- .rcppeigen.trans$ar2_unbounded(X)
    checkEquals(res$abs,  abs(X))
    checkEquals(res$abs2, X * X)
    checkEquals(res$exp,  exp(X))
    checkEquals(res$cos,  cos(X))
}

