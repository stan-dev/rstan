// load Rcpp
#include <Rcpp.h>
#include <Rmath.h>
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::Map;

//------------------------------------------------------------------------------
// unnormalized standard multivariate normal density function (log)
//------------------------------------------------------------------------------

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double log_densityRcpp_mu(NumericVector x, SEXP data, NumericVector mu) {

  VectorXd xe(as<Map<VectorXd> >(x));
  VectorXd mue(as<Map<VectorXd> >(mu));
  return -0.5*(xe - mue).transpose()*(xe - mue);

}
