#ifndef STAN_MATH_PRIM_SCAL_FUN_ROUND_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ROUND_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <boost/math/special_functions/round.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(double x) {
  if (is_nan(x))
    return std::numeric_limits<double>::quiet_NaN();
  return boost::math::round(x, boost_policy_t());
}

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(int x) { return round(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
