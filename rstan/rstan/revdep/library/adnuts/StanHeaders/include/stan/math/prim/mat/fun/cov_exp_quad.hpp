#ifndef STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of length scale
 *
 * @param x std::vector of elements that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
                 const T_l& length_scale) {
  using std::exp;
  check_positive("cov_exp_quad", "marginal variance", sigma);
  check_positive("cov_exp_quad", "length-scale", length_scale);
  for (size_t n = 0; n < x.size(); ++n)
    check_not_nan("cov_exp_quad", "x", x[n]);

  Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x.size(), x.size());

  int x_size = x.size();
  if (x_size == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);
  T_l neg_half_inv_l_sq = -0.5 / square(length_scale);

  for (int j = 0; j < (x_size - 1); ++j) {
    cov(j, j) = sigma_sq;
    for (int i = j + 1; i < x_size; ++i) {
      cov(i, j)
          = sigma_sq * exp(squared_distance(x[i], x[j]) * neg_half_inv_l_sq);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of std::vector of length scale
 *
 * @param x std::vector of elements that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma standard deviation
 * @param length_scale std::vector length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
                 const std::vector<T_l>& length_scale) {
  using std::exp;

  check_positive("cov_exp_quad", "marginal variance", sigma);
  for (size_t n = 0; n < x.size(); ++n) {
    check_not_nan("cov_exp_quad", "x", x[n]);
  }
  for (size_t n = 0; n < length_scale.size(); ++n) {
    check_positive("cov_exp_quad", "length-scale", length_scale[n]);
    check_not_nan("cov_exp_quad", "length-scale", length_scale[n]);
  }

  Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x.size(), x.size());

  int x_size = x.size();
  int l_size = length_scale.size();
  if (x_size == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);
  T_l temp_exp;

  for (int j = 0; j < (x_size - 1); ++j) {
    cov(j, j) = sigma_sq;
    for (int i = j + 1; i < x_size; ++i) {
      temp_exp = 0;
      for (int k = 0; k < l_size; ++k) {
        temp_exp += squared_distance(x[i], x[j]) / square(length_scale[k]);
      }
      cov(i, j) = sigma_sq * exp(-0.5 * temp_exp);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of of length scale
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const T_l& length_scale) {
  using std::exp;
  check_positive("cov_exp_quad", "marginal variance", sigma);
  check_positive("cov_exp_quad", "length-scale", length_scale);
  for (size_t n = 0; n < x1.size(); ++n)
    check_not_nan("cov_exp_quad", "x1", x1[n]);
  for (size_t n = 0; n < x2.size(); ++n)
    check_not_nan("cov_exp_quad", "x2", x2[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1.size(), x2.size());
  if (x1.size() == 0 || x2.size() == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);
  T_l neg_half_inv_l_sq = -0.5 / square(length_scale);

  for (size_t i = 0; i < x1.size(); ++i) {
    for (size_t j = 0; j < x2.size(); ++j) {
      cov(i, j)
          = sigma_sq * exp(squared_distance(x1[i], x2[j]) * neg_half_inv_l_sq);
    }
  }
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of length scale
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma standard deviation
 * @param length_scale std::vector of length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const std::vector<T_l>& length_scale) {
  using std::exp;
  check_positive("cov_exp_quad", "marginal variance", sigma);
  for (size_t n = 0; n < x1.size(); ++n)
    check_not_nan("cov_exp_quad", "x1", x1[n]);
  for (size_t n = 0; n < x2.size(); ++n)
    check_not_nan("cov_exp_quad", "x2", x2[n]);
  for (size_t n = 0; n < length_scale.size(); ++n)
    check_positive("cov_exp_quad", "length-scale", length_scale[n]);
  for (size_t n = 0; n < length_scale.size(); ++n)
    check_not_nan("cov_exp_quad", "length-scale", length_scale[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1.size(), x2.size());
  if (x1.size() == 0 || x2.size() == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);
  T_l temp_exp;

  for (size_t i = 0; i < x1.size(); ++i) {
    for (size_t j = 0; j < x2.size(); ++j) {
      temp_exp = 0;
      for (size_t k = 0; k < length_scale.size(); ++k) {
        temp_exp += squared_distance(x1[i], x2[j]) / square(length_scale[k]);
      }
      cov(i, j) = sigma_sq * exp(-0.5 * temp_exp);
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
