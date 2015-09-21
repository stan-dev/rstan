
#include <stan/math/prim/mat/fun/csr_extract_u.hpp>
#include <stan/math/prim/mat/fun/csr_extract_v.hpp>
#include <stan/math/prim/mat/fun/csr_extract_w.hpp>

#include <Rcpp.h>
#include <RcppEigen.h>

RcppExport SEXP extract_sparse_components(SEXP A) {
  Eigen::SparseMatrix<double> AA(Rcpp::as<Eigen::SparseMatrix<double> >(A));
  Eigen::SparseMatrix<double, Eigen::RowMajor> BB(AA.transpose());
  Eigen::Matrix<double, Eigen::Dynamic, 1> w = stan::math::csr_extract_w(BB);
  std::vector<double> ww(w.size());
  for (unsigned int i=0; i < w.size(); ++i) ww[i] = w.coeff(i);

  std::vector<int> v = stan::math::csr_extract_v(BB);
  std::vector<int> u = stan::math::csr_extract_u(BB);

  return Rcpp::List::create(
    Rcpp::Named("w") = ww,
    Rcpp::Named("v") = v,
    Rcpp::Named("u") = u
  );

}
