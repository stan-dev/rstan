#ifndef STAN_MATH_PRIM_MAT_PROB_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_POISSON_LOG_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Poisson distribution and log link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * @tparam T_n type of vector of variates (labels), integers >=0;
 * this can also be a single positive integer;
 * @tparam T_x type of the matrix of covariates (features); this
 * should be an Eigen::Matrix type whose number of rows should match the
 * length of n and whose number of columns should match the length of beta
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @tparam T_alpha type of the intercept;
 * this should be a single value;
 * @param n positive integer vector parameter
 * @param x design matrix
 * @param beta weight vector
 * @param alpha intercept (in log odds)
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if n is negative.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n, typename T_x, typename T_beta,
          typename T_alpha>
typename return_type<T_x, T_beta, T_alpha>::type poisson_log_glm_lpmf(
    const T_n &n, const T_x &x, const T_beta &beta, const T_alpha &alpha) {
  static const char *function = "poisson_log_glm_lpmf";
  typedef typename stan::partials_return_type<T_n, T_x, T_beta, T_alpha>::type
      T_partials_return;

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  if (!(stan::length(n) && stan::length(x) && stan::length(beta)))
    return 0.0;

  T_partials_return logp(0.0);

  check_nonnegative(function, "Vector of dependent variables", n);
  check_finite(function, "Matrix of independent variables", x);
  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_consistent_sizes(function, "Rows in matrix of independent variables",
                         x.col(0), "Vector of dependent variables", n);
  check_consistent_sizes(function, "Columns in matrix of independent variables",
                         x.row(0), "Weight vector", beta);

  if (!include_summand<propto, T_x, T_beta, T_alpha>::value)
    return 0.0;

  const size_t N = x.col(0).size();
  const size_t M = x.row(0).size();

  Matrix<T_partials_return, Dynamic, 1> n_vec(N, 1);
  {
    scalar_seq_view<T_n> n_seq_view(n);
    for (size_t n = 0; n < N; ++n) {
      n_vec[n] = n_seq_view[n];
    }
  }

  Matrix<T_partials_return, Dynamic, 1> beta_dbl(M, 1);
  {
    scalar_seq_view<T_beta> beta_vec(beta);
    for (size_t m = 0; m < M; ++m) {
      beta_dbl[m] = value_of(beta_vec[m]);
    }
  }
  Matrix<T_partials_return, Dynamic, Dynamic> x_dbl = value_of(x);

  Matrix<T_partials_return, Dynamic, 1> theta_dbl
      = (x_dbl * beta_dbl
         + Matrix<double, Dynamic, 1>::Ones(N, 1) * value_of(alpha));
  Matrix<T_partials_return, Dynamic, 1> exp_theta
      = theta_dbl.array().exp().matrix();

  for (size_t i = 0; i < N; i++) {
    // Compute the log-density.
    if (!(theta_dbl[i] == -std::numeric_limits<double>::infinity()
          && n_vec[i] == 0)) {
      if (include_summand<propto>::value)
        logp -= lgamma(n_vec[i] + 1.0);
      if (include_summand<propto, T_partials_return>::value)
        logp += n_vec[i] * theta_dbl[i] - exp_theta[i];
    }
  }

  // Compute the necessary derivatives.
  operands_and_partials<T_x, T_beta, T_alpha> ops_partials(x, beta, alpha);
  if (!(is_constant_struct<T_x>::value && is_constant_struct<T_beta>::value
        && is_constant_struct<T_alpha>::value)) {
    Matrix<T_partials_return, Dynamic, 1> theta_derivative = n_vec - exp_theta;
    if (!is_constant_struct<T_beta>::value) {
      ops_partials.edge2_.partials_ = x_dbl.transpose() * theta_derivative;
    }
    if (!is_constant_struct<T_x>::value) {
      ops_partials.edge1_.partials_ = theta_derivative * beta_dbl.transpose();
    }
    if (!is_constant_struct<T_alpha>::value) {
      ops_partials.edge3_.partials_[0] = theta_derivative.trace();
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_x, typename T_beta, typename T_alpha>
inline typename return_type<T_x, T_beta, T_alpha>::type poisson_log_glm_lpmf(
    const T_n &n, const T_x &x, const T_beta &beta, const T_alpha &alpha) {
  return poisson_log_glm_lpmf<false>(n, x, beta, alpha);
}
}  // namespace math
}  // namespace stan
#endif
