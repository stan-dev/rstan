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
typedef Eigen::PartialPivLU<Mat>        PPLU;
typedef Eigen::ColPivHouseholderQR<Mat> CPQR;
'

definitions <- list(
    "dense_PPLU" = list(signature(A_="matrix", b_="numeric"),
    '
    MMat           A(as<MMat>(A_));
    MVec           b(as<MVec>(b_));
    PPLU           lu(A);
    Mat            Ainv(lu.inverse());
    Vec            x(lu.solve(b));

    return List::create(Named("A",    A),
                        Named("Ainv", Ainv),
                        Named("b",    b),
                        Named("x",    x));
    '),
    "dense_CPQR" = list(signature(A_="matrix", b_="numeric"),
    '
    MMat           A(as<MMat>(A_));
    MVec           b(as<MVec>(b_));
    CPQR           qr(A);
    Mat            Ainv(qr.inverse());
    Vec            x(qr.solve(b));
    return List::create(Named("Ainv", Ainv),
                        Named("x",    x));
    ')
    )


.setUp <- function() {
    suppressMessages(require(inline))
    suppressMessages(require(RcppEigen))
    cxxargs <- ifelse(Rcpp:::capabilities()[["initializer lists"]],
                      "-std=c++0x","")
    tests <- ".rcppeigen.solve"
    if( ! exists( tests, globalenv() )) {
        fun <- RcppEigen:::compile_unit_tests(definitions,
                                              includes=incl,
                                              cxxargs = cxxargs)
        names(fun) <- names(definitions)
        assign( tests, fun, globalenv() )
    }
}

test.smallDense <- function() {
    A <- matrix(c(1,2,3,4), nrow=2L)
    B <- matrix(c(5,6,7,8), nrow=2L)
    b <- c(1,1)

    ## solutions to dense systems
    res <- .rcppeigen.solve$dense_PPLU(A, b)
    checkEquals(res$Ainv,  solve(A))
    checkEquals(res$x,     solve(A, b))

    res <- .rcppeigen.solve$dense_CPQR(A, b)
    checkEquals(res$Ainv,  solve(A))
    checkEquals(res$x,     solve(A, b))
}

test.largeDense <- function() {
    set.seed(1234321)
    N <- 100L
    AA <- matrix(rnorm(N * N), nrow=N)
    bb <- rnorm(N)

    res <- .rcppeigen.solve$dense_PPLU(AA, bb)
    checkEquals(res$Ainv,  solve(AA))
    checkEquals(res$x,     solve(AA, bb))

    res <- .rcppeigen.solve$dense_CPQR(AA, bb)
    checkEquals(res$Ainv,  solve(AA))
    checkEquals(res$x,     solve(AA, bb))
}

