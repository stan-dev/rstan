#ifndef STAN_MATH_PRIM_SCAL_FUN_TRUNC_HPP
#define STAN_MATH_PRIM_SCAL_FUN_TRUNC_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the nearest integral value that is not larger in
 * magnitude than the specified argument.
 *
 * @param[in] x Argument.
 * @return The truncated argument.
 */
inline double trunc(double x) {
  if (is_nan(x))
    return std::numeric_limits<double>::quiet_NaN();
  return boost::math::trunc(x, boost_policy_t());
}

/**
 * Return the nearest integral value that is not larger in
 * magnitude than the specified argument.
 *
 * @param[in] x Argument.
 * @return The truncated argument.
 */
inline double trunc(int x) { return trunc(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
