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
double log_densityRcpp(NumericVector x, SEXP data) {

  VectorXd xe(as<Map<VectorXd> >(x));
  return -0.5*xe.transpose()*xe;

}
