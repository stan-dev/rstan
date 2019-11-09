#ifndef STAN_MATH_PRIM_SCAL_FUN_ATANH_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ATANH_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/atanh.hpp>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic tangent of the specified value.
 * An argument of -1 returns negative infinity and an argument of 1
 * returns infinity.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is not in [-1, 1].
 */
inline double atanh(double x) {
  if (unlikely(is_nan(x)))
    return x;
  else
    return boost::math::atanh(x, boost_policy_t());
}

/**
 * Integer version of atanh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double atanh(int x) { return atanh(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
