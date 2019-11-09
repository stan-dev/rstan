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

test.wrap.R <- function(){

    fx <- cxxfunction( , '

    List vecs = List::create(
        _["Vec<complex>"] = Eigen::VectorXcd::Zero(5),
        _["Vec<double>"]  = Eigen::VectorXd::Zero(5),
        _["Vec<float>"]   = Eigen::VectorXf::Zero(5),
        _["Vec<int>"]     = Eigen::VectorXi::Zero(5)
    );

    // A VectorX<T> behaves as a matrix with one column but is converted to
    // a vector object in R, not a matrix of one column.  The distinction is
    // that VectorX<T> objects are defined at compile time to have one column,
    // whereas a MatrixX<T> has a dynamic number of columns that is set to 1
    // during execution of the code.  A MatrixX<T> object can be resized to have
    // a different number of columns.  A VectorX<T> object cannot.
    List cols = List::create(
        _["Col<complex>"] = Eigen::MatrixXcd::Zero(5, 1),
        _["Col<double>"]  = Eigen::MatrixXd::Zero(5, 1),
        _["Col<float>"]   = Eigen::MatrixXf::Zero(5, 1),
        _["Col<int>"]     = Eigen::MatrixXi::Zero(5, 1)
    );

    List rows = List::create(
        _["Row<complex>"] = Eigen::RowVectorXcd::Zero(5),
        _["Row<double>"]  = Eigen::RowVectorXd::Zero(5),
        _["Row<float>"]   = Eigen::RowVectorXf::Zero(5),
        _["Row<int>"]     = Eigen::RowVectorXi::Zero(5)
    );

    List matrices = List::create(
        _["Mat<complex>"] = Eigen::MatrixXcd::Identity(3, 3),
        _["Mat<double>"]  = Eigen::MatrixXd::Identity(3, 3),
        _["Mat<float>"]   = Eigen::MatrixXf::Identity(3, 3),
        _["Mat<int>"]     = Eigen::MatrixXi::Identity(3, 3)
    );

    // ArrayXX<t> objects have the same structure as matrices but allow
    // componentwise arithmetic.  A * B is matrix multiplication for
    // matrices and componentwise multiplication for arrays.
    List arrays2 = List::create(
        _["Arr2<complex>"] = Eigen::ArrayXXcd::Zero(3, 3),
        _["Arr2<double>"]  = Eigen::ArrayXXd::Zero(3, 3),
        _["Arr2<float>"]   = Eigen::ArrayXXf::Zero(3, 3),
        _["Arr2<int>"]     = Eigen::ArrayXXi::Zero(3, 3)
    );

    // ArrayX<t> objects have the same structure as VectorX<T> objects
    // but allow componentwise arithmetic, including functions like exp, log,
    // sqrt, ...
    List arrays1 = List::create(
        _["Arr1<complex>"] = Eigen::ArrayXcd::Zero(5),
        _["Arr1<double>"]  = Eigen::ArrayXd::Zero(5),
        _["Arr1<float>"]   = Eigen::ArrayXf::Zero(5),
        _["Arr1<int>"]     = Eigen::ArrayXi::Zero(5)
    );

    List operations = List::create(
        _["Op_seq"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10),  // arguments are length.out, start, end
        _["Op_log"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).log(),
        _["Op_exp"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).exp(),
        _["Op_sqrt"] = Eigen::ArrayXd::LinSpaced(6, 1, 10).sqrt(),
        _["Op_cos"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).cos()
    );

    List output = List::create(
    	_["vectors : VectorX<T>"]   = vecs,
    	_["matrices : MatrixX<T>"]  = matrices,
    	_["rows : RowVectorX<T>"]   = rows,
    	_["columns : MatrixX<T>"]   = cols,
        _["arrays2d : ArrayXX<T>"]  = arrays2,
        _["arrays1d : ArrayX<T>"]   = arrays1,
        _["operations : ArrayXd"]   = operations
        );

    return output ;
	' , plugin = "RcppEigen" )

    res <- fx()

    checkEquals( res[["vectors : VectorX<T>"]][["Vec<complex>"]], complex(5), msg = "VectorXcd::Zero(5)")
    checkEquals( res[["vectors : VectorX<T>"]][["Vec<double>"]], double(5), msg = "VectorXd::Zero(5)")
    checkEquals( res[["vectors : VectorX<T>"]][["Vec<float>"]], double(5), msg = "VectorXf::Zero(5)")
    checkEquals( res[["vectors : VectorX<T>"]][["Vec<int>"]], integer(5), msg = "VectorXi::Zero(5)")

    checkEquals( res[["matrices : MatrixX<T>"]][["Mat<complex>"]], (1+0i) * diag(nr=3L), msg = "MatrixXcd::Identity(3,3)")
    checkEquals( res[["matrices : MatrixX<T>"]][["Mat<double>"]], diag(nr=3L), msg = "MatrixXd::Identity(3,3)")
    checkEquals( res[["matrices : MatrixX<T>"]][["Mat<float>"]], diag(nr=3L), msg = "MatrixXf::Identity(3,3)")
    checkEquals( res[["matrices : MatrixX<T>"]][["Mat<int>"]], matrix(as.integer((diag(nr=3L))),nr=3L), msg = "MatrixXi::Identity(3,3)")

    checkEquals( res[["rows : RowVectorX<T>"]][["Row<complex>"]], matrix(complex(5), nr=1L), msg = "RowVectorXcd::Zero(5)")
    checkEquals( res[["rows : RowVectorX<T>"]][["Row<double>"]], matrix(numeric(5), nr=1L), msg = "RowVectorXd::Zero(5)")
    checkEquals( res[["rows : RowVectorX<T>"]][["Row<float>"]], matrix(numeric(5), nr=1L), msg = "RowVectorXf::Zero(5)")
    checkEquals( res[["rows : RowVectorX<T>"]][["Row<int>"]], matrix(integer(5), nr=1L), msg = "RowVectorXi::Zero(5)")

    checkEquals( res[["columns : MatrixX<T>"]][["Col<complex>"]], as.matrix(complex(5)), msg = "MatrixXcd::Zero(5, 1)")
    checkEquals( res[["columns : MatrixX<T>"]][["Col<double>"]], as.matrix(numeric(5)), msg = "MatrixXd::Zero(5, 1)")
    checkEquals( res[["columns : MatrixX<T>"]][["Col<float>"]], as.matrix(numeric(5)), msg = "MatrixXf::Zero(5, 1)")
    checkEquals( res[["columns : MatrixX<T>"]][["Col<int>"]], as.matrix(integer(5)), msg = "MatrixXi::Zero(5, 1)")

    checkEquals( res[["arrays2d : ArrayXX<T>"]][["Arr2<complex>"]], matrix(complex(9L), nc=3L), msg = "ArrayXXcd::Zero(3,3)")
    checkEquals( res[["arrays2d : ArrayXX<T>"]][["Arr2<double>"]], matrix(numeric(9L), nc=3L), msg = "ArrayXXd::Zero(3,3)")
    checkEquals( res[["arrays2d : ArrayXX<T>"]][["Arr2<float>"]], matrix(numeric(9L), nc=3L), msg = "ArrayXXf::Zero(3,3)")
    checkEquals( res[["arrays2d : ArrayXX<T>"]][["Arr2<int>"]], matrix(integer(9L), nc=3L), msg = "ArrayXXi::Zero(3,3)")

    checkEquals( res[["arrays1d : ArrayX<T>"]][["Arr1<complex>"]], complex(5), msg = "ArrayXcd::Zero(5)")
    checkEquals( res[["arrays1d : ArrayX<T>"]][["Arr1<double>"]], double(5), msg = "ArrayXd::Zero(5)")
    checkEquals( res[["arrays1d : ArrayX<T>"]][["Arr1<float>"]], double(5), msg = "ArrayXf::Zero(5)")
    checkEquals( res[["arrays1d : ArrayX<T>"]][["Arr1<int>"]], integer(5), msg = "ArrayXi::Zero(5)")

    oneTen <- seq(1, 10, length.out=6L)

    checkEquals( res[["operations : ArrayXd"]][["Op_seq"]],  oneTen,       msg = "Op_seq")
    checkEquals( res[["operations : ArrayXd"]][["Op_log"]],  log(oneTen),  msg = "Op_log")
    checkEquals( res[["operations : ArrayXd"]][["Op_exp"]],  exp(oneTen),  msg = "Op_exp")
    checkEquals( res[["operations : ArrayXd"]][["Op_sqrt"]], sqrt(oneTen), msg = "Op_sqrt")
    checkEquals( res[["operations : ArrayXd"]][["Op_cos"]],  cos(oneTen),  msg = "Op_cos")

}

test.as.Vec <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    Eigen::VectorXi                                m1 = input[0] ; /* implicit as */
    Eigen::VectorXd                                m2 = input[1] ; /* implicit as */
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> m3 = input[0] ; /* implicit as */
    Eigen::VectorXf                                m4 = input[1] ; /* implicit as */

    List res = List::create(m1.sum(), m2.sum(), m3.sum(), m4.sum());

    return res ;

    ', plugin = "RcppEigen" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Vec>" )
}

test.as.MVec <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    const Eigen::Map<Eigen::VectorXi>   m1 = input[0] ; // maps share storage and do not allow conversion
    const Eigen::Map<Eigen::VectorXd>   m2 = input[1] ;

    List res = List::create(m1.sum(), m2.sum());

    return res ;

    ', plugin = "RcppEigen" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 2 ), msg = "as<MVec>" )
}

test.as.MRowVec <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    const Eigen::Map<Eigen::RowVectorXi>   m1 = input[0] ; // maps share storage and do not allow conversion
    const Eigen::Map<Eigen::RowVectorXd>   m2 = input[1] ;

    List res = List::create(m1.sum(), m2.sum());

    return res ;

    ', plugin = "RcppEigen" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 2 ), msg = "as<MRowVec>" )
}


test.as.MMat <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    const Eigen::Map<Eigen::MatrixXi>   m1 = input[0]; // maps share storage and do not allow conversion
    const Eigen::Map<Eigen::MatrixXd>   m2 = input[1] ;
// FIXME: Write a version of as specifically for complex matrices.
//    const Eigen::Map<Eigen::MatrixXcd>  m3 = input[2] ;

    List res = List::create(m1.sum(), m2.sum());//, m3.sum());

    return res ;

    ', plugin = "RcppEigen" )

    integer_mat <- matrix(as.integer(diag(nr=4L)), nc=4L)
    numeric_mat <- diag(nr=5L)
    complex_mat <- (1+0i) * diag(nr=5L)
    res <- fx(list(integer_mat, numeric_mat, complex_mat))
    checkEquals(unlist(res), c(4L, 5)#, 5+0i)
                , msg = "as<MMat>" )
}

test.as.MSpMat <- function() {
    suppressMessages(require("Matrix"))
    data("KNex", package = "Matrix")
    fx <- cxxfunction( signature(input_ = "list"), '
    List input(input_) ;
    const Eigen::MappedSparseMatrix<double>  m1 = input[0]; // maps share storage and do not allow conversion

    List res = List::create(_["nnz"]   = double(m1.nonZeros()),
                            _["nr"]    = double(m1.rows()),
                            _["nc"]    = double(m1.cols()),
                            _["inSz"]  = double(m1.innerSize()),
                            _["outSz"] = double(m1.outerSize()),
                            _["sum"]   = m1.sum());

    return res ;

    ', plugin = "RcppEigen" )

    KNX <- KNex[[1]]
    res <- fx(KNex)
    checkEquals(unname(unlist(res)),
                as.numeric(c(nnzero(KNX), nrow(KNX), ncol(KNX), nrow(KNX), ncol(KNX), sum(KNX@x))),
                msg = "as<MSPMatrix>")
}

test.as.SpMat <- function() {
    suppressMessages(require("Matrix"))
    data("KNex", package = "Matrix")
    fx <- cxxfunction( signature(input_ = "list"), '
    List input(input_) ;
    const Eigen::SparseMatrix<double>  m1 = input[0];

    List res = List::create(_["nnz"]   = double(m1.nonZeros()),
                            _["nr"]    = double(m1.rows()),
                            _["nc"]    = double(m1.cols()),
                            _["inSz"]  = double(m1.innerSize()),
                            _["outSz"] = double(m1.outerSize()),
                            _["sum"]   = m1.sum());

    return res ;
    ', plugin = "RcppEigen" )

    KNX <- KNex[[1]]
    res <- fx(KNex)
    checkEquals(unname(unlist(res)),
                as.numeric(c(nnzero(KNX), nrow(KNX), ncol(KNX), nrow(KNX), ncol(KNX), sum(KNX@x))),
                msg = "as<MSPMatrix>")
}
