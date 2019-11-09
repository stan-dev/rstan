#!/usr/bin/r -t
#
# Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
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
    suppressMessages(require(inline))
}

test.wrapSparse.double.R <- function(){

    fx <- cxxfunction( , '

    Eigen::SparseMatrix<double>  mm(9,3);
    mm.reserve(9);
    for (int j = 0; j < 3; ++j) {
        mm.startVec(j);
        for (int i = 3 * j; i < 3 * (j + 1); ++i)
            mm.insertBack(i, j) = 1.;
    }
    mm.finalize();
    return wrap(mm);
' , plugin = "RcppEigen" )

    res <- fx()
    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    checkEquals( res, rr, msg = "wrap<SparseMatrix<double> >")
}

test.wrapSparse.double.ColMajor.R <- function(){

    fx <- cxxfunction( , '

    Eigen::SparseMatrix<double, Eigen::ColMajor>  mm(9,3);
    mm.reserve(9);
    for (int j = 0; j < 3; ++j) {
        mm.startVec(j);
        for (int i = 3 * j; i < 3 * (j + 1); ++i)
            mm.insertBack(i, j) = 1.;
    }
    mm.finalize();
    return wrap(mm);
' , plugin = "RcppEigen" )

    res <- fx()
    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    checkEquals( res, rr, msg = "wrap<SparseMatrix<double, Eigen::ColMajor> >")
}

## test.wrapSparse.int.ColMajor.R <- function(){  ## classes not yet exported from Matrix

##     fx <- cxxfunction( , '

##     Eigen::SparseMatrix<int, Eigen::ColMajor>  mm(9,3);
##     mm.reserve(9);
##     for (int j = 0; j < 3; ++j) {
##         mm.startVec(j);
##         for (int i = 3 * j; i < 3 * (j + 1); ++i)
##             mm.insertBack(i, j) = 1;
##     }
##     mm.finalize();
##     return wrap(mm);
## ' , plugin = "RcppEigen" )

##     #res <- fx()
##     #rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
##     #colnames(rr) <- NULL
##     #checkEquals( res, rr, msg = "wrap<SparseMatrix<double, Eigen::ColMajor> >")
##     checkException( fx(), msg = "wrap<SparseMatrix<int, Eigen::ColMajor> >" )
## }

test.wrapSparse.double.RowMajor.R <- function(){

    fx <- cxxfunction( , '

    Eigen::SparseMatrix<double, Eigen::RowMajor>  mm(9,3);
    mm.reserve(9);
    for (int irow = 0; irow < 9; ++irow) {
        mm.startVec(irow);
        mm.insertBack(irow, irow / 3) = static_cast<double>( 9 - irow );
    }
    mm.finalize();
    return wrap(mm);
' , plugin = "RcppEigen" )

    res <- fx()
    rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
    colnames(rr) <- NULL
    checkEquals( res, rr, msg = "wrap<SparseMatrix<double, Eigen::RowMajor> >")
}

## test.wrapSparse.int.RowMajor.R <- function(){

##     fx <- cxxfunction( , '

##     Eigen::SparseMatrix<int, Eigen::RowMajor>  mm(9,3);
##     mm.reserve(9);
##     for (int irow = 0; irow < 9; ++irow) {
##         mm.startVec(irow);
##         mm.insertBack(irow, irow / 3) = 9 - irow;
##     }
##     mm.finalize();
##     return wrap(mm);
## ' , plugin = "RcppEigen" )

##     #res <- fx()
##     #rr <- new( "igRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=9L:1L, Dim=c(9L,3L) )
##     #colnames(rr) <- NULL
##     #checkEquals( res, rr, msg = "wrap<SparseMatrix<int, Eigen::RowMajor> >")
##     checkException( fx(), msg = "wrap<SparseMatrix<int, Eigen::RowMajor> >" )
## }

test.asSparse.double.ColMajor.R <- function(){

    fx <- cxxfunction( sig=signature(R_mm="dgCMatrix"), '

    Eigen::SparseMatrix<double, Eigen::ColMajor> mm = Rcpp::as<Eigen::SparseMatrix<double, Eigen::ColMajor> >( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<SparseMatrix<double, Eigen::ColMajor> >")
}

test.asMappedSparse.double.ColMajor.R <- function(){

    fx <- cxxfunction( sig=signature(R_mm="dgCMatrix"), '

    typedef Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor> > MapMat;
    MapMat mm = Rcpp::as<MapMat>( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<Map<SparseMatrix<double, Eigen::ColMajor> > >")
}

test.asMappedSparse.deprecated.double.ColMajor.R <- function(){

    fx <- cxxfunction( sig=signature(R_mm="dgCMatrix"), '
    // Deprecated
    typedef Eigen::MappedSparseMatrix<double, Eigen::ColMajor> MapMat;
    MapMat mm = Rcpp::as<MapMat>( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<MappedSparseMatrix<double, Eigen::ColMajor> >")
}

test.asSparse.double.RowMajor.R <- function(){
    fx <- cxxfunction( sig=signature(R_mm="dgRMatrix"), '

    Eigen::SparseMatrix<double, Eigen::RowMajor> mm = Rcpp::as<Eigen::SparseMatrix<double, Eigen::RowMajor> >( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<SparseMatrix<double, Eigen::RowMajor> >")
}

test.asMappedSparse.double.RowMajor.R <- function(){
    fx <- cxxfunction( sig=signature(R_mm="dgRMatrix"), '

    typedef Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor> > MapMat;
    MapMat mm = Rcpp::as<MapMat>( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<Map<SparseMatrix<double, Eigen::RowMajor> > >")
}

test.asMappedSparse.deprecated.double.RowMajor.R <- function(){
    fx <- cxxfunction( sig=signature(R_mm="dgRMatrix"), '
    // Deprecated
    typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MapMat;
    MapMat mm = Rcpp::as<MapMat>( R_mm );
    return wrap(mm);
' , plugin = "RcppEigen" )

    rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
    colnames(rr) <- NULL
    res <- fx( R_mm = rr )
    checkEquals( res, rr, msg = "as<MappedSparseMatrix<double, Eigen::RowMajor> >")
}


test.sparseCholesky.R <- function() {
    suppressMessages(require("Matrix", character.only=TRUE))
    data("KNex", package = "Matrix")

    fx <- cxxfunction( signature(input_ = "list"), '
    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::SparseMatrix;
    using Eigen::SimplicialLDLT;
    using Eigen::Success;

    List input(input_);
    const Map<SparseMatrix<double> > m1 = input[0];
    const Map<VectorXd>              v1 = input[1];
    SparseMatrix<double>             m2(m1.cols(), m1.cols());
    m2.selfadjointView<Lower>().rankUpdate(m1.adjoint());

    SimplicialLDLT<SparseMatrix<double> > ff(m2);
    VectorXd                        res = ff.solve(m1.adjoint() * v1);

    return List::create(_["res"]   = res,
                        _["rows"]  = double(ff.rows()),
                        _["cols"]  = double(ff.cols()));
',
                      plugin = "RcppEigen")

    rr <- fx(KNex)
    checkEquals(rr[[1]], as.vector(solve(crossprod(KNex[[1]]),
                                         crossprod(KNex[[1]], KNex[[2]])),
                                   mode="numeric"),
                "Cholmod solution")
}
